"""
An iterator owning an open Stockholm stream and, for URLs, its temporary download.
Construct it with [`eachmsa`](@ref). Iteration consumes the stream and cannot restart it.
"""
mutable struct StockholmMSAIterator{T,S<:IO,K}
    io::S
    kwargs::K
    temporary::Union{Nothing,String}
    closed::Bool
    ready::Bool
end

Base.IteratorSize(::Type{<:StockholmMSAIterator}) = Base.SizeUnknown()
Base.eltype(::Type{<:StockholmMSAIterator{T}}) where {T} = T
Base.isopen(msas::StockholmMSAIterator) = !msas.closed && isopen(msas.io)

function Base.close(msas::StockholmMSAIterator)
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

# Look ahead only as far as the next header. In particular, isempty and zip must
# not consume an alignment from this stateful iterator.
function Base.isdone(msas::StockholmMSAIterator, ::Nothing = nothing)
    msas.closed && return true
    msas.ready && return false
    try
        while !eof(msas.io)
            line = strip(readline(msas.io))
            isempty(line) && continue
            line == "# STOCKHOLM 1.0" ||
                throw(ArgumentError("Expected a # STOCKHOLM 1.0 header, got: $line"))
            msas.ready = true
            return false
        end
        close(msas)
        return true
    catch
        close(msas)
        rethrow()
    end
end

function Base.iterate(msas::StockholmMSAIterator{T}, ::Nothing = nothing) where {T}
    Base.isdone(msas) && return nothing
    try
        msas.ready = false
        msa = parse_file(msas.io, Stockholm, T; msas.kwargs...)
        return msa, nothing
    catch
        close(msas)
        rethrow()
    end
end

"""
    eachmsa(source, Stockholm[, output::Type]; kwargs...)
    eachmsa(f, source, Stockholm[, output::Type]; kwargs...)

Iterate over the alignments in a Stockholm file, parsing one MSA at a time. The default
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
    ::Type{Stockholm},
    ::Type{T} = AnnotatedMultipleSequenceAlignment;
    kwargs...,
) where {T}
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
        msas = StockholmMSAIterator{T,typeof(io),typeof(options)}(
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
    format::Type{Stockholm},
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

# TODO: Support concatenated Clustal alignments, using each new CLUSTAL header as
# the next MSA boundary while keeping wrapped blocks within the same alignment.
