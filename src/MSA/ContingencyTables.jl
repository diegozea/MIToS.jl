type ContingencyTable{T,N,A} <: AbstractArray{T,N}
    alphabet::ResidueAlphabet
    temporal::Array{T,N}
    table::NamedArray{T,N}
    marginals::NamedArray{T,2}
    total::T
end
