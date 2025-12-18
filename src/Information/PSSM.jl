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

@inline function _collect_background(background, alphabet::ResidueAlphabet)
    values =
        if background isa Probabilities || background isa ContingencyTable
            gettablearray(background)
        else
            background
        end
    vector = Float64.(vec(values))
    length(vector) == length(alphabet) ||
        throw(ArgumentError("Background length $(length(vector)) doesn't match alphabet length $(length(alphabet))."))
    total = sum(vector)
    total > 0 || throw(ArgumentError("Background distribution sums to zero."))
    if total != 1
        vector ./= total
    end
    vector
end

@inline function _log_odds(
    p::Float64,
    q::Float64,
    zero_policy::Symbol,
    use_log2::Bool,
    invlogbase::Float64,
    epsT::Float64,
)
    if zero_policy === :clamp
        p = max(p, epsT)
        q = max(q, epsT)
    elseif zero_policy === :error
        (p > 0 && q > 0) || throw(DomainError((p, q), "Zero probability encountered."))
    else
        # zero_policy === :negInf
        if p == 0 && q == 0
            return NaN
        elseif p == 0
            return -Inf
        elseif q == 0
            return Inf
        end
    end
    ratio = p / q
    use_log2 ? log2(ratio) : log(ratio) * invlogbase
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
- `background = BLOSUM62_Pi`: background distribution `q(a)`; accepts vectors,
  `Probabilities` or `ContingencyTable` objects. It is normalized if needed.
- `base::Number = 2`: logarithm base (e.g., 2 for bits).
- `gap_handling::Symbol = :skip`: if `:skip`, gaps are ignored when the alphabet includes
  them; if `:include`, gaps are scored using their observed frequency.
- `zero_policy::Symbol = :negInf`: handling of zero probabilities. `:negInf` produces
  `-Inf`, `Inf` or `NaN`; `:error` throws; `:clamp` replaces zeros with `eps`.
- `all_gap_value = NaN`: value used for every residue in columns with no valid
  observations after applying the gap policy.

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
    gap_handling::Symbol = :skip,
    zero_policy::Symbol = :negInf,
    all_gap_value = NaN,
)
    dims == 2 || throw(ArgumentError("pssm supports dims=2 (per-column) only."))
    gap_handling in (:skip, :include) ||
        throw(ArgumentError("Unsupported gap_handling = $(gap_handling)."))
    zero_policy in (:negInf, :error, :clamp) ||
        throw(ArgumentError("Unsupported zero_policy = $(zero_policy)."))
    base_val = Float64(base)
    (base_val > 0) && (base_val != 1) ||
        throw(ArgumentError("Base must be positive and different from 1."))

    q = _collect_background(background, alphabet)

    nres = length(alphabet)
    ncols = size(msa, 2)
    scores = Matrix{Float64}(undef, nres, ncols)

    column_table = ContingencyTable(Float64, Val{1}, alphabet)
    gap_index = alphabet[GAP]
    include_gap = gap_handling === :include || gap_index > nres
    use_log2 = base_val == 2
    invlogbase = use_log2 ? 1.0 : inv(log(base_val))
    epsT = eps(Float64)

    freq_array = gettablearray(column_table)

    @inbounds for j = 1:ncols
        cleanup!(column_table)
        col_view = @view msa[:, j]
        frequencies!(column_table, col_view; weights = weights, pseudocounts = pseudocounts)

        if !include_gap && gap_index <= nres
            freq_array[gap_index] = 0.0
            update_marginals!(column_table)
        end

        total = gettotal(column_table)
        if total == 0
            scores[:, j] .= all_gap_value
            continue
        end

        normalize!(column_table)
        prob_view = gettablearray(column_table)
        for i = 1:nres
            p = prob_view[i]
            qi = q[i]
            scores[i, j] = _log_odds(p, qi, zero_policy, use_log2, invlogbase, epsT)
        end
    end

    PSSMResult(scores, alphabet, q, base_val)
end
