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
    Utils._close_input(msas.io, msas.temporary)
end

"""
Read one prepared alignment, forwarding the output type and parsing keywords. Formats
whose parser consumes the next header can specialize this method and set `msas.ready = true`
to preserve that header for the next iteration.
"""
function _read_msa!(msas::MSAIterator{F,T}) where {F,T}
    parse_file(msas.io, F, T; msas.kwargs...)
end

function _read_msa!(msas::MSAIterator{Clustal,T}) where {T}
    IDS, SEQS, annot, has_next =
        _load_clustal_sequences(eachline(msas.io); header_read = true)
    msas.ready = has_next
    _parse_msa((IDS, SEQS, annot), T; msas.kwargs...)
end

# Look ahead only as far as the next header. In particular, isempty and zip must
# not consume an alignment from this stateful iterator.
function Base.isdone(msas::MSAIterator{F}, ::Nothing = nothing) where {F}
    msas.closed && return true
    msas.ready && return false
    try
        msas.ready = _read_msa_header(msas.io, _msa_header(F); strict = true)
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
    eachmsa(source, format::Type[, output::Type]; kwargs...)
    eachmsa(f, source, format::Type[, output::Type]; kwargs...)

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

[`read_file`](@ref) reads only the first alignment and warns if another alignment is found.
"""
function eachmsa(
    source::AbstractString,
    ::Type{F},
    ::Type{T} = AnnotatedMultipleSequenceAlignment;
    kwargs...,
) where {F<:Union{Stockholm,Clustal},T}
    remote = Utils._is_url(source)
    temporary = remote ? Utils._download_tempname(source) : nothing
    filename = temporary === nothing ? source : temporary
    io = nothing
    try
        if remote
            download_file(source, filename; headers = Dict("Accept-Encoding" => "identity"))
        end
        io = open(filename, "r")
        io = Utils._input_stream(io, source)
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
        Utils._close_input(io, temporary)
        rethrow()
    end
end

function eachmsa(
    f::Function,
    source::AbstractString,
    format::Type{<:Union{Stockholm,Clustal}},
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

"""
Read the first MSA and warn if another alignment header is found. If the parser leaves
`has_next` unchecked (`nothing`), look for the next header after parsing.
"""
function Utils._read_file(
    io::IO,
    format::Type{F},
    output::Type{T} = AnnotatedMultipleSequenceAlignment;
    kwargs...,
) where {F<:Union{Stockholm,Clustal},T}
    IDS, SEQS, annot, has_next = _load_sequences(io, format, output)
    msa = _parse_msa((IDS, SEQS, annot), output; kwargs...)
    if has_next === nothing
        has_next = _read_msa_header(io, _msa_header(format))
    end
    if has_next
        @warn "Read only the first alignment; use `eachmsa` to read all."
    end
    msa
end
