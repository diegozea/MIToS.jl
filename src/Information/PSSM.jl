# Position-Specific Scoring Matrix (PSSM)
# =======================================

"""
    struct PSSMResult{T,A}
        scores::Matrix{T}
        alphabet::A
        background::Vector{Float64}
        base::Float64
    end

Result returned by [`pssm`](@ref), containing the log-odds scores for a protein PSSM.
`scores` is a matrix whose rows follow the provided `alphabet` ordering and whose columns
match the alignment positions.
"""
struct PSSMResult{T,A}
    scores::Matrix{T}
    alphabet::A
    background::Vector{Float64}
    base::Float64
end

function _collect_background(background::AbstractArray, alphabet::ResidueAlphabet)
    values = _gettablearray(background)
    n_values = length(values)
    n_ab = length(alphabet)
    if n_values != n_ab
        throw(
            ArgumentError(
                "Background length $n_values doesn't match alphabet length $n_ab.",
            ),
        )
    end
    vector = Vector{Float64}(undef, n_ab)
    total = 0.0
    @inbounds for (i, val) in enumerate(values)
        value = Float64(val)
        isfinite(value) ||
            throw(DomainError(value, "Background values must be finite Float64 values."))
        value ≥ 0 || throw(
            DomainError(value, "Background values must be nonnegative Float64 values."),
        )
        vector[i] = value
        total += value
    end
    if !(total > 0)
        throw(DomainError(total, "Background distribution must have positive sum."))
    end
    if total != 1.0
        vector ./= total
    end
    vector
end

"""
    pssm(msa::AbstractArray{Residue}; kwargs...) -> PSSMResult

Compute a log-odds position-specific scoring matrix (PSSM) from a protein MSA using
MIToS probability estimation. Scores are returned as a matrix with rows ordered by the
provided alphabet and columns corresponding to alignment positions.

# Keywords

  - `dims::Int = 2`: compute scores per column. Other values throw an `ArgumentError`.
  - `alphabet = UngappedAlphabet()`: residue alphabet; score rows follow this order.
  - `weights = NoClustering()`: sequence weights used during probability estimation.
  - `pseudocounts = NoPseudocount()`: pseudocounts applied before normalization.
  - `background = BLOSUM62_Pi`: background distribution `q(a)`; accepts `AbstractArray`,
    `Probabilities` or `ContingencyTable` objects. It is normalized if needed.
  - `base::Number = 2`: logarithm base (e.g., 2 for bits).

Gaps are handled by the chosen `alphabet`: `UngappedAlphabet()` ignores gaps, while
`GappedAlphabet()` includes them.

Scores are computed directly as log-odds, so IEEE floating-point rules apply for zeros
(`-Inf`, `Inf`, or `NaN`).

Columns with no valid observations (e.g. all gaps using `UngappedAlphabet()`) are filled
with `NaN`.

# Examples

```julia
using MIToS.Information, MIToS.MSA

msa = hcat(res"AC", res"AD", res"AE")' # three sequences, two positions
result = pssm(msa; background = fill(1 / 20, 20))
```
"""
function pssm(
    msa::AbstractArray{Residue};
    dims::Int = 2,
    alphabet::ResidueAlphabet = UngappedAlphabet(),
    weights::WeightTypes = NoClustering(),
    pseudocounts::Pseudocount = NoPseudocount(),
    background = BLOSUM62_Pi,
    base::Number = 2,
)
    dims == 2 || throw(ArgumentError("pssm supports dims=2 (per-column) only."))
    base_val = Float64(base)
    (base_val > 0) && (base_val != 1) ||
        throw(ArgumentError("Base must be positive and different from 1."))

    q = _collect_background(background, alphabet)

    nres = length(alphabet)
    ncols = size(msa, 2)
    scores = Matrix{Float64}(undef, nres, ncols)

    column_table = ContingencyTable(Float64, Val{1}, alphabet)
    use_log2 = base_val == 2
    invlogbase = use_log2 ? 1.0 : inv(log(base_val))

    @inbounds for j = 1:ncols
        cleanup!(column_table)
        col_view = @view msa[:, j]
        frequencies!(column_table, col_view; weights = weights, pseudocounts = pseudocounts)

        total = gettotal(column_table)
        if total == 0
            scores[:, j] .= NaN
            continue
        end

        normalize!(column_table)
        prob_view = gettablearray(column_table)
        for i = 1:nres
            p = prob_view[i]
            qi = q[i]
            ratio = p / qi
            scores[i, j] = use_log2 ? log2(ratio) : log(ratio) * invlogbase
        end
    end

    PSSMResult(scores, alphabet, q, base_val)
end
