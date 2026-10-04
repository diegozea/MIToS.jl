"""
Validate collections handled by the generic writer, leaving custom writers in control
of their own inputs.
"""
function Utils._validate_write(
    msas,
    format::Type{F},
    ::Type{S},
) where {F<:Union{Stockholm,Clustal},S<:IO}
    msas isa AbstractMatrix{Residue} && return nothing
    which(print_file, (S, typeof(msas), Type{F})) ===
    which(print_file, (S, Any, Type{F})) || return nothing
    Base.IteratorEltype(typeof(msas)) isa Base.HasEltype ||
        throw(ArgumentError("Iterator must declare an alignment eltype."))
    T = eltype(msas)
    T <: AbstractMatrix{Residue} ||
        throw(ArgumentError("Expected an alignment eltype; got $T."))
    nothing
end

"""
    print_file(io::IO, msas, format::Type; kwargs...)

Write a collection or iterator of alignments in `Stockholm` or `Clustal` format.
The iterator must declare an element type compatible with `AbstractMatrix{Residue}`.
Alignments are written one at a time, passing any format keywords to each alignment's
writer. An empty typed collection writes nothing. The caller owns the output stream and
any input stream.
"""
function Utils.print_file(io::IO, msas, format::Type{<:Union{Stockholm,Clustal}}; kwargs...)
    Utils._validate_write(msas, format, typeof(io))
    for msa in msas
        print_file(io, msa::AbstractMatrix{Residue}, format; kwargs...)
    end
    nothing
end
