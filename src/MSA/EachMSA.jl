function hasnextmsa(
    io::TranscodingStream,
    format::Type{<:Union{Stockholm,Clustal}};
    strict::Bool = false,
)
    hasnextmsa(io, _msa_header(format); strict = strict)
end

support_eachmsa(::Type{Stockholm}) = true
support_eachmsa(::Type{Clustal}) = true

"""
Reuse a buffered input, or add a buffer to an ordinary input.
"""
_buffer_msa_input(io::IO) = NoopStream(io)
_buffer_msa_input(io::TranscodingStream) = io

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

# Look ahead only as far as the next header. In particular, isempty and zip must
# not consume an alignment from this stateful iterator.
function Base.isdone(msas::MSAIterator{F}, ::Nothing = nothing) where {F}
    msas.closed && return true
    msas.ready && return false
    try
        msas.ready = hasnextmsa(msas.io, F; strict = true)
        msas.ready && return false
        close(msas)
        return true
    catch
        close(msas)
        rethrow()
    end
end

function Base.iterate(msas::MSAIterator{F,T}, ::Nothing = nothing) where {F,T}
    Base.isdone(msas) && return nothing
    try
        msas.ready = false
        msa = parse_file(msas.io, F, T; msas.kwargs...)
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
Each alignment is read with `parse_file`, including user-defined output types.
Additional formats can declare [`support_eachmsa`](@ref) as `true` and implement
[`hasnextmsa`](@ref) and `parse_file` without defining another iterator type.

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

eachmsa(\"Pfam-A.seed.gz\", Stockholm) do msas
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
) where {F<:MSAFormat,T}
    support_eachmsa(F) || throw(MethodError(eachmsa, (source, F, T)))
    remote = Utils._is_url(source)
    temporary = remote ? Utils._download_tempname(source) : nothing
    filename = temporary === nothing ? source : temporary
    io = nothing
    try
        if remote
            download_file(source, filename; headers = Dict("Accept-Encoding" => "identity"))
        end
        io = open(filename, "r")
        io = _buffer_msa_input(Utils._input_stream(io, source))
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

"""
Read one alignment through `parse_file` and warn if another alignment header follows.
For other formats, forward the original arguments to preserve their parser's defaults.
"""
function Utils._read_file(
    io::IO,
    format::Type{F},
    args::Vararg{Any,N};
    kwargs...,
) where {F<:MSAFormat,N}
    support_eachmsa(format) || return parse_file(io, format, args...; kwargs...)
    io = _buffer_msa_input(io)
    output_args = isempty(args) ? (AnnotatedMultipleSequenceAlignment,) : args
    msa = parse_file(io, format, output_args...; kwargs...)
    if hasnextmsa(io, format)
        @warn "Read only the first alignment; use `eachmsa` to read all."
    end
    msa
end
