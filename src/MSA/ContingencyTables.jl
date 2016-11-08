type ContingencyTable{T,N,A} <: AbstractArray{T,N}
    alphabet::ResidueAlphabet
    temporal::Array{T,N}
    table::NamedArray{T,N}
    marginals::NamedArray{T,2}
    total::T
end

# AbstractArray Interface
# -----------------------

Base.size(table::ContingencyTable) = size(table.table)

Base.getindex(table::ContingencyTable, i...) = getindex(table.table, i...)

function Base.getindex(t::ContingencyTable, i::Residue)
    index = t.alphabet[i]
    @assert index != 22 "Residue $i isn't in the alphabet"
    getindex(t.table, index)
end

function Base.getindex(table::ContingencyTable, i::Residue, j::Residue)
    index_i = t.alphabet[i]
    index_j = t.alphabet[j]
    @assert index_i != 22 "Residue $i isn't in the alphabet"
    @assert index_j != 22 "Residue $j isn't in the alphabet"
    getindex(table.table, index_i, index_j)
end

function Base.getindex(t::ContingencyTable, I::Residue...)
    N = length(I)
    indexes = Array(Int,N)
    for i in 1:N
        index = t.alphabet[I[i]]
        @assert index != 22 "Residue $i isn't in the alphabet"
        indexes[i] = index
    end
    getindex(t.table, indexes...)
end
