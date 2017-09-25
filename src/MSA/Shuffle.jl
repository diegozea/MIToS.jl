"""
It's like `Base.shuffle`. When a `Matrix{Residue}` is used, you can indicate if the gaps
should remain their positions using the last boolean argument. The previous argument should
be the dimension to shuffle, 1 for shuffling residues in a sequence (row) or 2 for shuffling
residues in a column.

```julia
julia> msa = hcat(res"RRE",res"DDK", res"G--")
3×3 Array{MIToS.MSA.Residue,2}:
 R  D  G
 R  D  -
 E  K  -

julia> srand(42);

julia> shuffle(msa, 1, true)
3×3 Array{MIToS.MSA.Residue,2}:
 G  D  R
 D  R  -
 E  K  -

julia> srand(42);

julia> shuffle(msa, 1, false)
3×3 Array{MIToS.MSA.Residue,2}:
 G  D  R
 D  -  R
 -  E  K

```
"""
function Base.shuffle!(r::AbstractRNG, msa::Matrix{Residue},
                       dim::Int, fixedgaps::Bool=true)
    nseq, ncol = size(msa)
    @assert dim == 1 || dim == 2 "The dimension must be 1 (sequences) or 2 (columns)"
    if fixedgaps
        mask = msa .!= GAP
        if dim == 2
            @inbounds for i in 1:ncol
                shuffle!(view(msa,mask[:,i],i))
            end
        elseif dim == 1
            @inbounds for i in 1:nseq
                shuffle!(view(msa,i,mask[i,:]))
            end
        end
    else
        if dim == 2
            @inbounds for i in 1:ncol
                shuffle!(view(msa,:,i))
            end
        elseif dim == 1
            @inbounds for i in 1:nseq
                shuffle!(view(msa,i,:))
            end
        end
    end
    msa
end

function Base.shuffle!(msa::Matrix{Residue}, args...)
    shuffle!(GLOBAL_RNG, msa, args...)
end

"""
It's like `shuffle` but in-place. When a `Matrix{Residue}` or a `AbstractAlignedObject`
(sequence or MSA) is used, you can indicate if the gaps should remain their positions
using the last boolean argument.
"""
function Base.shuffle(r::AbstractRNG, msa::Matrix{Residue}, args...)
    shuffle!(r, copy(msa), args...)
end

function Base.shuffle(msa::Matrix{Residue}, args...)
    shuffle!(GLOBAL_RNG, copy(msa), args...)
end

function Base.shuffle(r::AbstractRNG,
                      msa::Union{AbstractAlignedObject, NamedArray{Residue,2}}, args...)
    shuffle(r, copy(getresidues(msa)), args...)
end

function Base.shuffle(msa::Union{AbstractAlignedObject, NamedArray{Residue,2}}, args...)
    shuffle(GLOBAL_RNG, copy(getresidues(msa)), args...)
end
