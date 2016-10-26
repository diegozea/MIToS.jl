"""
When a `Matrix{Residue}` is used, you can indicate if the gaps should remain their
positions using the last boolean argument.
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
When a `Matrix{Residue}` or a `AbstractAlignedObject` (sequence or MSA) is used, you can
indicate if the gaps should remain their positions using the last boolean argument.
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
