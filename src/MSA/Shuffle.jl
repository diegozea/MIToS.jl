function _subset_indices(msa::Matrix{Residue}, dims::Int, subset, fixed_reference)
    _indices = if subset === Colon()
        nseq, ncol = size(msa)
        if dims == 1
            1:nseq
        else
            1:ncol
        end
    else
        if eltype(subset) !== Int
            throw(
                ArgumentError(
                    "For a Matrix{Residue}, subset must be an iterator of Int values or Colon()",
                ),
            )
        end
        if isa(subset, AbstractRange)
            collect(subset)
        elseif isa(subset, Int)
            [subset]
        else
            subset
        end
    end
    indices = convert(Vector{Int}, _indices)
    if fixed_reference && dims == 1
        filter!(!=(1), indices)
    end
    indices
end

function _subset_indices(msa::NamedResidueMatrix, dims::Int, subset, fixed_reference)
    dict = dims == 1 ? msa.dicts[1] : msa.dicts[2]
    _indices = NamedArrays.indices(dict, subset)
    indices = convert(Vector{Int}, isa(_indices, Int) ? [_indices] : _indices)
    if fixed_reference && dims == 1
        filter!(!=(1), indices)
    end
    indices
end

function _subset_indices(
    msa::AbstractMultipleSequenceAlignment,
    dims,
    subset,
    fixed_reference,
)::Vector{Int}
    _subset_indices(msa.matrix, dims, subset, fixed_reference)
end

function shuffle_msa!(
    r::AbstractRNG,
    msa::AbstractMatrix{Residue},
    subset = Colon();
    dims::Int = 2,
    fixedgaps::Bool = true,
    fixed_reference::Bool = false,
)
    @assert dims == 1 || dims == 2 "dims must be 1 for shuffling along sequences or 2 for columns"
    subset_indices = _subset_indices(msa, dims, subset, fixed_reference)
    msa_matrix = getresidues(msa)
    nseq, ncol = size(msa_matrix)
    mask = fixedgaps ? msa_matrix .!= GAP : trues(nseq, ncol)
    if fixed_reference
        mask[1, :] .= 0
    end
    for i in subset_indices
        to_shuffle =
            dims == 1 ? view(msa_matrix, i, mask[i, :]) : view(msa_matrix, mask[:, i], i)
        shuffle!(r, to_shuffle)
    end
    msa
end


function shuffle_msa!(r::AbstractRNG, msa::MultipleSequenceAlignment, args...; kwargs...)
    shuffle_msa!(r, msa.matrix, args...; kwargs...)
    msa
end

function shuffle_msa!(
    r::AbstractRNG,
    msa::AnnotatedMultipleSequenceAlignment,
    subset = Colon();
    dims::Int = 2,
    fixedgaps::Bool = true,
    fixed_reference::Bool = false,
)
    shuffle_msa!(r, msa.matrix, subset; dims, fixedgaps, fixed_reference)

    # Annotate the modifications
    subset_indices = _subset_indices(msa, dims, subset, fixed_reference)
    n = length(subset_indices)
    entities = dims == 1 ? "sequences" : "columns"
    message = "$n $entities shuffled."
    fixed = if fixedgaps && fixed_reference
        " Gaps and residues in the first sequence"
    elseif fixedgaps
        " Gaps"
    elseif fixed_reference
        " Residues in the first sequence"
    else
        ""
    end
    if !isempty(fixed)
        message *= fixed
        message *= " were kept in their positions."
    end
    annotate_modification!(msa.annotations, message)
    if dims == 1
        seqnames = sequencenames(msa)
        for i in subset_indices
            seqname = seqnames[i]
            setannotsequence!(msa, seqname, "Shuffled", "true")
            # Delete SeqMap of the shuffled sequences
            delete!(msa.annotations.sequences, (seqname, "SeqMap"))
        end
    else
        shuffled = zeros(Int, ncolumns(msa))
        shuffled[subset_indices] .= 1
        setannotcolumn!(msa, "Shuffled", join(shuffled))
    end
    msa
end


shuffle_msa_doc = md"""
It randomly permute residues in the MSA `msa` along sequences (`dims=1`) or columns 
(`dims=2`, the default). The optional positional argument `subset` allows to shuffle only 
a subset of the sequences or columns. The optional keyword argument `fixedgaps` indicates 
if the gaps should remain their positions (`true` by default). The optional keyword
argument `fixed_reference` indicates if the residues in the first sequence should remain
in their positions (`false` by default).
"""

"""
    shuffle_msa!([rng=default_rng(),] msa::AbstractMatrix{Residue}, subset=Colon(); dims=2, fixedgaps=true, fixed_reference=false)

In-place version of [`shuffle_msa`](@ref). $shuffle_msa_doc
"""
function shuffle_msa!(msa::AbstractMatrix{Residue}, args...; kwargs...)
    shuffle_msa!(Random.default_rng(), msa, args...; kwargs...)
end

function shuffle_msa(r::AbstractRNG, msa::AbstractMatrix{Residue}, args...; kwargs...)
    shuffle_msa!(r, deepcopy(msa), args...; kwargs...)
end

"""
    shuffle_msa([rng=default_rng(),] msa::AbstractMatrix{Residue}, subset=Colon(); dims=2, fixedgaps=true, fixed_reference=false)

$shuffle_msa_doc To shuffle in-place, see [`shuffle_msa!`](@ref).

```jldoctest
julia> using MIToS.MSA

julia> using Random

julia> msa = hcat(res"RRE",res"DDK", res"G--")
3×3 Matrix{Residue}:
 R  D  G
 R  D  -
 E  K  -

julia> Random.seed!(42);

julia> shuffle_msa(msa, dims=1, fixedgaps=true)
3×3 Matrix{Residue}:
 G  D  R
 R  D  -
 E  K  -

julia> Random.seed!(42);

julia> shuffle_msa(msa, dims=1, fixedgaps=false)
3×3 Matrix{Residue}:
 G  D  R
 R  -  D
 E  K  -

```
"""
function shuffle_msa(msa::AbstractMatrix{Residue}, args...; kwargs...)
    shuffle_msa(Random.default_rng(), msa, args...; kwargs...)
end

"""
It's like `Random.shuffle`. When a `Matrix{Residue}` is used, you can indicate if the gaps
should remain their positions using the last boolean argument. The previous argument should
be the dimension to shuffle, 1 for shuffling residues in a sequence (row) or 2 for shuffling
residues in a column.

**DEPRECATED:** This method is deprecated. Use [`shuffle_msa!`](@ref) instead.
"""
function Random.shuffle!(
    r::AbstractRNG,
    msa::AbstractMatrix{Residue},
    dim::Int,
    fixedgaps::Bool = true,
)
    @warn "The function `shuffle!(r, msa, dim, fixedgaps)` is deprecated. Use `shuffle_msa!(r, msa; dims, fixedgaps)` instead."
    shuffle_msa!(r, msa, Colon(); dims = dim, fixedgaps = fixedgaps) |> getresidues
end

function Random.shuffle!(msa::AbstractMatrix{Residue}, args...)
    shuffle!(Random.default_rng(), msa, args...)
end

"""
It's like `shuffle` but in-place. When a `Matrix{Residue}` or a `AbstractAlignedObject`
(sequence or MSA) is used, you can indicate if the gaps should remain their positions
using the last boolean argument.

**DEPRECATED:** This method is deprecated. Use [`shuffle_msa`](@ref) instead.
"""
function Random.shuffle(
    r::AbstractRNG,
    msa::AbstractMatrix{Residue},
    dim::Int,
    fixedgaps::Bool = true,
)
    @warn "The function `shuffle(r, msa, dim, fixedgaps)` is deprecated. Use `shuffle_msa(r, msa; dims, fixedgaps)` instead."
    shuffle_msa(r, msa, Colon(); dims = dim, fixedgaps = fixedgaps) |> getresidues
end

function Random.shuffle(msa::AbstractMatrix{Residue}, args...)
    shuffle(Random.GLOBAL_RNG, msa, args...)
end
