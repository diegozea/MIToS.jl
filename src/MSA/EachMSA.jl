"""
An iterator owning an open alignment stream and, for URLs, its temporary download.
Construct it with [`eachmsa`](@ref). Iteration consumes the stream and cannot restart it.
The format parameter selects how to find and parse the next alignment.
"""
mutable struct MSAIterator{F<:MSAFormat,T,S<:IO,K}
    io::S
    kwargs::K
    temporary::Union{Nothing,String}
    closed::Bool
    ready::Bool
end

Base.IteratorSize(::Type{<:MSAIterator}) = Base.SizeUnknown()
Base.eltype(::Type{<:MSAIterator{F,T}}) where {F,T} = T
Base.isopen(msas::MSAIterator) = !msas.closed && isopen(msas.io)

function Base.close(msas::MSAIterator)
    msas.closed && return nothing
    msas.closed = true
    try
        close(msas.io)
    finally
        if msas.temporary !== nothing
            rm(msas.temporary; force = true)
        end
    end
    nothing
end

"""
Check format support before opening a stream or downloading a URL. Formats supporting
iteration define this method and `_prepare_msa!`, and provide a parser for one alignment.
"""
function _check_eachmsa_format(::Type{F}) where {F<:MSAFormat}
    throw(ArgumentError("eachmsa does not support the $F format"))
end

_check_eachmsa_format(::Type{Stockholm}) = nothing
_check_eachmsa_format(::Type{Clustal}) = nothing

"""
Prepare the next alignment without parsing it. Return `true` if an alignment is ready
to parse, or `false` at the end of the stream. The iterator caches this result so repeated
lookahead does not consume an alignment.
"""
function _prepare_msa!(msas::MSAIterator{Stockholm})
    while !eof(msas.io)
        line = strip(readline(msas.io))
        isempty(line) && continue
        line == "# STOCKHOLM 1.0" ||
            throw(ArgumentError("Expected a # STOCKHOLM 1.0 header, got: $line"))
        return true
    end
    false
end

function _prepare_msa!(msas::MSAIterator{Clustal})
    while !eof(msas.io)
        line = strip(readline(msas.io))
        isempty(line) && continue
        _is_clustal_header(line) ||
            throw(ArgumentError("Expected a CLUSTAL header, got: $line"))
        return true
    end
    false
end

"""
Read one prepared alignment, forwarding the output type and parsing keywords. Formats
whose `parse_file` consumes more than one alignment can specialize this method. If reading
an alignment also prepares the next one, set `msas.ready = true` to preserve that lookahead.
"""
function _read_msa!(msas::MSAIterator{F,T}) where {F,T}
    parse_file(msas.io, F, T; msas.kwargs...)
end

function _read_msa!(msas::MSAIterator{Clustal,T}) where {T}
    # A new header ends the current alignment; blank lines only separate its blocks.
    # Save that lookahead for the next iteration, including on non-seekable gzip streams.
    lines = Iterators.takewhile(eachline(msas.io)) do line
        if _is_clustal_header(line)
            msas.ready = true
            return false
        end
        true
    end
    _parse_msa(T; msas.kwargs...) do create_annotations
        _load_clustal_sequences(lines)
    end
end

# Look ahead only as far as the next header. In particular, isempty and zip must
# not consume an alignment from this stateful iterator.
function Base.isdone(msas::MSAIterator, ::Nothing = nothing)
    msas.closed && return true
    msas.ready && return false
    try
        msas.ready = _prepare_msa!(msas)
        msas.ready && return false
        close(msas)
        return true
    catch
        close(msas)
        rethrow()
    end
end

function Base.iterate(msas::MSAIterator, ::Nothing = nothing)
    Base.isdone(msas) && return nothing
    try
        msas.ready = false
        msa = _read_msa!(msas)
        return msa, nothing
    catch
        close(msas)
        rethrow()
    end
end

"""
    eachmsa(source, format[, output::Type]; kwargs...)
    eachmsa(f, source, format[, output::Type]; kwargs...)

Iterate over the alignments in a `Stockholm` or `Clustal` file, parsing one MSA at a time.
For Clustal, each alignment starts with its own CLUSTAL header. The default
output is `AnnotatedMultipleSequenceAlignment`; the output types and parsing keywords
are the same as for [`read_file`](@ref) and [`parse_file`](@ref).

The source is a local path or an HTTP, HTTPS or FTP URL. Files ending in `.gz` are
decompressed incrementally using one open stream. A URL is downloaded once to a temporary
file when `eachmsa` is called. No alignment is parsed until iteration begins, and previous
alignments are not retained by the iterator. Each individual MSA must still fit in memory.

The iterator is consumed as it is read: iterating again continues from its current position.
An empty file yields no MSAs, and a file with one alignment yields once.

Exhaustion or a parsing error closes the stream and removes any temporary download.
Use a `do` block to ensure cleanup when your analysis finishes, even if you stop early
or your analysis throws an exception. It returns the result of the block.

```julia
using MIToS.MSA

eachmsa(\"Pfam-A.full.gz\", Stockholm) do msas
    for msa in msas
        println(getannotfile(msa, \"AC\", \"\"), '\\t', nsequences(msa))
    end
end
```

[`read_file`](@ref) returns the first alignment of a Stockholm file as a single MSA object.
"""
function eachmsa(
    source::AbstractString,
    ::Type{F},
    ::Type{T} = AnnotatedMultipleSequenceAlignment;
    kwargs...,
) where {F<:MSAFormat,T}
    _check_eachmsa_format(F)
    remote = any(prefix -> startswith(source, prefix), ("http://", "https://", "ftp://"))
    temporary = remote ? tempname() * (endswith(source, ".gz") ? ".gz" : "") : nothing
    filename = temporary === nothing ? source : temporary
    io = nothing
    try
        if remote
            download_file(source, filename; headers = Dict("Accept-Encoding" => "identity"))
        end
        io = open(filename, "r")
        if endswith(source, ".gz")
            io = GzipDecompressorStream(io)
        end
        options = (; kwargs...)
        msas = MSAIterator{F,T,typeof(io),typeof(options)}(
            io,
            options,
            temporary,
            false,
            false,
        )
        finalizer(close, msas)
        return msas
    catch
        try
            io === nothing || close(io)
        finally
            temporary === nothing || rm(temporary; force = true)
        end
        rethrow()
    end
end

function eachmsa(
    f::Function,
    source::AbstractString,
    format::Type{<:MSAFormat},
    args...;
    kwargs...,
)
    msas = eachmsa(source, format, args...; kwargs...)
    try
        return f(msas)
    finally
        close(msas)
    end
end
