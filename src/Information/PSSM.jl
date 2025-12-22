# Position-Specific Scoring Matrix (PSSM)
# =======================================

# Abstract Types
# ==============

abstract type AbstractPositionSpecificMatrix{T,N,A} <: AbstractArray{T,N} end

# Result Type
# ===========

"""
    struct PositionSpecificScoreMatrix{T,A}
        table::NamedArray{T,2,Array{T,2},NTuple{2,OrderedDict{String,Int}}}
        alphabet::A
        base::Union{Irrational{:ℯ}, Float64, Int}
    end

Result returned by [`position_specific_scoring_matrix`](@ref), containing the log-odds scores for a protein PSSM.
`table` is a `NamedArray` whose rows follow the provided `alphabet` ordering and whose columns
match the alignment positions.
"""
struct PositionSpecificScoreMatrix{T,A} <: AbstractPositionSpecificMatrix{T,2,A}
    table::NamedArray{T,2,Array{T,2},NTuple{2,OrderedDict{String,Int}}}
    alphabet::A
    base::Union{Irrational{:ℯ}, Float64, Int}
end

# Getters
# -------

@inline getalphabet(scores::AbstractPositionSpecificMatrix) = scores.alphabet
@inline gettable(scores::AbstractPositionSpecificMatrix) = scores.table
@inline gettablearray(scores::AbstractPositionSpecificMatrix) = getarray(gettable(scores))

# AbstractArray
# -------------

for f in (:size, :getindex)
    @eval Base.$(f)(scores::AbstractPositionSpecificMatrix, args...) =
        $(f)(gettable(scores), args...)
end

# Show
# ----

function Base.show(
    io::IO,
    ::MIME"text/plain",
    scores::PositionSpecificScoreMatrix,
)
    base_val = scores.base
    unit =
        base_val == ℯ ? "nats" : base_val == 2 ? "bits" : base_val == 10 ? "hartleys" :
        nothing
    print(io, typeof(scores), " log base=", base_val)
    if unit !== nothing
        print(io, " units=", unit)
    end
    print(io, " : ")
    print(io, "\ntable : ")
    show(io, MIME"text/plain"(), gettable(scores))
end

# Base Conversion
# ---------------

@inline _convert_base(base::Irrational{:ℯ}) = base
@inline _convert_base(base::Int) = base

@inline function _convert_base(base::Integer)
    if base >= 0 && base <= typemax(Int)
        return Int(base)
    end
    Float64(base)
end

@inline _convert_base(base::Number) = Float64(base)

# Background Collection
# ---------------------

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

# PSSM
# ====

"""
    position_specific_scoring_matrix(msa::AbstractArray{Residue}; kwargs...) -> PositionSpecificScoreMatrix

Compute a log-odds position-specific scoring matrix (PSSM) from a protein MSA using
MIToS probability estimation. Scores are returned as a `NamedArray` with rows ordered by
the provided alphabet and columns corresponding to alignment positions.

# Keywords

  - `alphabet = UngappedAlphabet()`: residue alphabet; score rows follow this order.
  - `weights = NoClustering()`: sequence weights used during probability estimation.
  - `pseudocounts = NoPseudocount()`: pseudocounts applied before normalization.
  - `background = BLOSUM62_Pi`: background distribution `q(a)`; accepts `AbstractArray`,
    `Probabilities` or `ContingencyTable` objects. It is normalized if needed.
  - `base::Number = ℯ`: logarithm base.

Use the keyword argument `base` to change the base of the log. $_DOC_LOG_BASE

Gaps are handled by the chosen `alphabet`: `UngappedAlphabet()` ignores gaps, while
`GappedAlphabet()` includes them.

Scores are computed directly as log-odds, so IEEE floating-point rules apply for zeros
(`-Inf`, `Inf`, or `NaN`).

Columns with no valid observations (e.g. all gaps using `UngappedAlphabet()`) are filled
with `NaN`.

# Examples

```julia
using MIToS.Information, MIToS.MSA

msa = permutedims(hcat(res"AC", res"AD", res"AE")) # three sequences, two positions
result = position_specific_scoring_matrix(msa; background = fill(1 / 20, 20))
```
"""
function position_specific_scoring_matrix(
    msa::AbstractArray{Residue};
    alphabet::ResidueAlphabet = UngappedAlphabet(),
    weights::WeightTypes = NoClustering(),
    pseudocounts::Pseudocount = NoPseudocount(),
    background = BLOSUM62_Pi,
    base::Number = ℯ,
)
    if !(base > 0) || base == 1 # this also catches NaN
        throw(
            ArgumentError(
                "The logarithm base must be positive and different from 1 (base=$base).",
            ),
        )
    end

    q = _collect_background(background, alphabet)

    nres = length(alphabet)
    ncols = ncolumns(msa)
    scores = Matrix{Float64}(undef, nres, ncols)

    column_table = ContingencyTable(Float64, Val{1}, alphabet)
    invlogbase = base === ℯ ? 1.0 : inv(log(base))

    for j = 1:ncols
        cleanup!(column_table)
        col_view = @inbounds @view msa[:, j]
        # counts:
        frequencies!(column_table, col_view; weights = weights, pseudocounts = pseudocounts)
        # probabilities:
        normalize!(column_table)
        P = gettablearray(column_table)
        # log-odds scores:
        for i = 1:nres
            p = @inbounds P[i]
            qi = @inbounds q[i]
            ratio = p / qi
            @inbounds scores[i, j] =
                invlogbase == 1.0 ? log(ratio) : log(ratio) * invlogbase
        end
    end

    row_dict = getnamedict(alphabet)
    col_names = columnnames(msa)
    col_dict = OrderedDict{String,Int}(col_names[i] => i for i = 1:ncols)
    table = NamedArray(scores, (row_dict, col_dict), ("Res", "Col"))

    PositionSpecificScoreMatrix(table, alphabet, _convert_base(base))
end
