"""
Accept a single alignment or a collection with a declared alignment element type.
"""
function Utils._validate_write(msas, format::Type{<:Union{Stockholm,Clustal}})
    msas isa AbstractMatrix{Residue} && return nothing
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
    Utils._validate_write(msas, format)
    for msa in msas
        print_file(io, msa::AbstractMatrix{Residue}, format; kwargs...)
    end
    nothing
end
