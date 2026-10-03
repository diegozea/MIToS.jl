# Position-Specific Scoring Matrix (PSSM)
# =======================================

# Docstring reusable parts
# ------------------------

const _LOG_BASE_UNITS = string(
    "Units follow the log base: base=2 yields bits, ",
    "base=ℯ yields nats, base=10 yields hartleys.",
)

const _DOC_PROFILE_SHAPE = string(
    "rows following the provided alphabet ordering and ",
    "columns corresponding to alignment positions.",
)

const _DOC_MSA_KWARGS = """
  - `alphabet = UngappedAlphabet()`: residue alphabet defining the row order.
  - `weights = NoClustering()`: sequence weights to reduce redundancy.
  - `pseudocounts = NoPseudocount()`: smoothing applied before normalization.
"""
const _DOC_PSSM_FORMULA = string(
    "``Sᵢ(a) = log_b(pᵢ(a) / q(a))``, where i is the alignment column, ",
    "a is a residue symbol in the alphabet, pᵢ(a) is the PSPM probability at ",
    "that column, q(a) is the background probability for residue a, and b is ",
    "the log base.",
)
const _DOC_PSPM_FORMULA = string(
    "``pᵢ(a) = fᵢ(a) / ∑ₐ fᵢ(a)``, where i is the alignment column, ",
    "a is a residue symbol in the alphabet, and ``∑ₐ`` sums over all residues in ",
    "that alphabet.",
)

# Abstract Types
# ==============

"""
    AbstractPositionSpecificMatrix{A}

Abstract supertype for position-wise residue profiles derived from a protein MSA. These
profiles are matrix-like with shape (n_residues × n_positions), where rows follow the
chosen residue alphabet and columns correspond to alignment positions.

Concrete subtypes represent three stages of a profile:

  - PFM ([`PositionFrequencyMatrix`](@ref)): ``fᵢ(a)`` is the (possibly weighted and
    pseudocounted) count of residue a in column i.
  - PSPM ([`PositionSpecificProbabilityMatrix`](@ref)): $(_DOC_PSPM_FORMULA)
  - PSSM ([`PositionSpecificScoringMatrix`](@ref)): $(_DOC_PSSM_FORMULA)

These profiles summarize column composition and can be used to score aligned sequences.
"""
abstract type AbstractPositionSpecificMatrix{A} <: AbstractArray{Float64,2} end

# Result Types
# ============

"""
    struct PositionFrequencyMatrix{A}
        table::NamedArray{Float64,2,Array{Float64,2},NTuple{2,OrderedDict{String,Int}}}
        alphabet::A
    end

Position frequency matrix (PFM) for a protein MSA representing the frequencies/counts of
each residue at each alignment position. So that ``fᵢ(a)`` is the (possibly weighted and
pseudocounted) count of residue a in column i. Weights can down-weight redundant sequences,
and pseudocounts smooth columns with sparse observations.
A `PositionFrequencyMatrix` is a subtype of [`AbstractPositionSpecificMatrix`](@ref). That
stores the frequency matrix in the `table` field. `table` is a `NamedArray` with
$(_DOC_PROFILE_SHAPE).
"""
mutable struct PositionFrequencyMatrix{A} <: AbstractPositionSpecificMatrix{A}
    table::NamedArray{Float64,2,Array{Float64,2},NTuple{2,OrderedDict{String,Int}}}
    alphabet::A
end

"""
    struct PositionSpecificProbabilityMatrix{A}
        table::NamedArray{Float64,2,Array{Float64,2},NTuple{2,OrderedDict{String,Int}}}
        alphabet::A
    end

Position-specific probability matrix (PSPM). Each column is a probability distribution
over the alphabet at a position: $(_DOC_PSPM_FORMULA) Columns are normalized independently;
if a column has zero total weight (for example, an all-gap column under an ungapped
alphabet), that column is filled with `NaN`.
"""
mutable struct PositionSpecificProbabilityMatrix{A} <: AbstractPositionSpecificMatrix{A}
    table::NamedArray{Float64,2,Array{Float64,2},NTuple{2,OrderedDict{String,Int}}}
    alphabet::A
end

"""
    struct PositionSpecificScoringMatrix{A}
        table::NamedArray{Float64,2,Array{Float64,2},NTuple{2,OrderedDict{String,Int}}}
        alphabet::A
        base::Union{Irrational{:ℯ}, Float64, Int}
    end

Position-specific scoring matrix (PSSM) with log-odds scores relative to a background
distribution. Scores are defined as $(_DOC_PSSM_FORMULA)

Positive scores indicate enrichment over background; negative scores indicate depletion.
The `base` field controls units. $(_LOG_BASE_UNITS) IEEE floating-point rules apply for
zeros (pᵢ(a) = 0 ⇒ `-Inf`, q(a) = 0 ⇒ `Inf`, invalid ratios ⇒ `NaN`).

To score an aligned sequence x of length n_positions, you can sum per-position scores:
``score(x) = ∑ᵢ Sᵢ(xᵢ)``, where xᵢ is the residue observed at position i.

`table` is a `NamedArray` with $(_DOC_PROFILE_SHAPE).
"""
mutable struct PositionSpecificScoringMatrix{A} <: AbstractPositionSpecificMatrix{A}
    table::NamedArray{Float64,2,Array{Float64,2},NTuple{2,OrderedDict{String,Int}}}
    alphabet::A
    base::Union{Irrational{:ℯ},Float64,Int}
end

# Constructors
# ============

# Helpers: Arguments
# ------------------

@inline function _check_matrix_nrows(data::AbstractMatrix, alphabet::ResidueAlphabet)
    nrows = size(data, 1)
    nrows == length(alphabet) || throw(
        ArgumentError(
            "Row count ($nrows) doesn't match alphabet length ($(length(alphabet))).",
        ),
    )
end

# Helpers: Log Base
# -----------------

@inline function _check_log_base(base::Number)
    if !(base > 0) || base == 1 # this also catches NaN
        throw(
            ArgumentError(
                "The logarithm base must be positive and different from 1 (base=$base).",
            ),
        )
    end
end

@inline _convert_base(base::Irrational{:ℯ}) = base
@inline _convert_base(base::Int) = base

@inline function _convert_base(base::Integer)
    if base >= 0 && base <= typemax(Int)
        return Int(base)
    end
    Float64(base)
end

@inline _convert_base(base::Number) = Float64(base)

# NamedArray constructors (main)
# ------------------------------

function PositionFrequencyMatrix(
    table::NamedArray{Float64,2,Array{Float64,2},NTuple{2,OrderedDict{String,Int}}},
    alphabet::A,
) where {A<:ResidueAlphabet}
    _check_matrix_nrows(table, alphabet)
    PositionFrequencyMatrix{A}(table, alphabet)
end

function PositionSpecificProbabilityMatrix(
    table::NamedArray{Float64,2,Array{Float64,2},NTuple{2,OrderedDict{String,Int}}},
    alphabet::A,
) where {A<:ResidueAlphabet}
    _check_matrix_nrows(table, alphabet)
    PositionSpecificProbabilityMatrix{A}(table, alphabet)
end

function PositionSpecificScoringMatrix(
    table::NamedArray{Float64,2,Array{Float64,2},NTuple{2,OrderedDict{String,Int}}},
    alphabet::A,
    base::Union{Irrational{:ℯ},Float64,Int},
) where {A<:ResidueAlphabet}
    _check_log_base(base)
    _check_matrix_nrows(table, alphabet)
    PositionSpecificScoringMatrix{A}(table, alphabet, base)
end

function PositionSpecificScoringMatrix(
    table::NamedArray{Float64,2,Array{Float64,2},NTuple{2,OrderedDict{String,Int}}},
    alphabet::A,
    base::Number = ℯ,
) where {A<:ResidueAlphabet}
    PositionSpecificScoringMatrix(table, alphabet, _convert_base(base))
end


# Matrix constructors (AbstractMatrix -> NamedArray)
# --------------------------------------------------

function _named_matrix(data::AbstractMatrix, alphabet::ResidueAlphabet)
    ncols = size(data, 2)
    row_dict = getnamedict(alphabet)
    col_dict = OrderedDict{String,Int}(string(i) => i for i = 1:ncols)
    NamedArray(Matrix{Float64}(data), (row_dict, col_dict), ("Res", "Col"))
end

PositionFrequencyMatrix(data::AbstractMatrix, alphabet::ResidueAlphabet) =
    PositionFrequencyMatrix(_named_matrix(data, alphabet), alphabet)

PositionSpecificProbabilityMatrix(data::AbstractMatrix, alphabet::ResidueAlphabet) =
    PositionSpecificProbabilityMatrix(_named_matrix(data, alphabet), alphabet)

function PositionSpecificScoringMatrix(
    data::AbstractMatrix,
    alphabet::ResidueAlphabet,
    base::Number = ℯ,
)
    PositionSpecificScoringMatrix(_named_matrix(data, alphabet), alphabet, base)
end

# Size-based constructors (ncol -> Matrix)
# ----------------------------------------

PositionFrequencyMatrix(alphabet::ResidueAlphabet, ncols::Int) =
    PositionFrequencyMatrix(zeros(Float64, length(alphabet), ncols), alphabet)

PositionSpecificProbabilityMatrix(alphabet::ResidueAlphabet, ncols::Int) =
    PositionSpecificProbabilityMatrix(fill(NaN, length(alphabet), ncols), alphabet)

function PositionSpecificScoringMatrix(
    alphabet::ResidueAlphabet,
    ncols::Int,
    base::Number = ℯ,
)
    PositionSpecificScoringMatrix(zeros(Float64, length(alphabet), ncols), alphabet, base)
end

# Convenience constructors
# ------------------------

Base.zeros(::Type{PositionFrequencyMatrix}, alphabet::ResidueAlphabet, ncols::Int) =
    PositionFrequencyMatrix(alphabet, ncols)

Base.zeros(
    ::Type{PositionSpecificScoringMatrix},
    alphabet::ResidueAlphabet,
    ncols::Int,
    base::Number = ℯ,
) = PositionSpecificScoringMatrix(alphabet, ncols, base)

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

Base.setindex!(scores::AbstractPositionSpecificMatrix, v, args...) =
    setindex!(gettable(scores), v, args...)

@inline function Base.getindex(scores::AbstractPositionSpecificMatrix, r::Residue, c::Int)
    i = scores.alphabet[r]
    table = gettablearray(scores)
    @boundscheck checkbounds(table, i, c)
    @inbounds table[i, c]
end

@inline function Base.setindex!(
    scores::AbstractPositionSpecificMatrix,
    v,
    r::Residue,
    c::Int,
)
    i = scores.alphabet[r]
    table = gettablearray(scores)
    @boundscheck checkbounds(table, i, c)
    @inbounds table[i, c] = v
end

# Show
# ----

function Base.show(io::IO, ::MIME"text/plain", scores::PositionSpecificScoringMatrix)
    base_val = scores.base
    unit =
        base_val == ℯ ? "nats" :
        base_val == 2 ? "bits" : base_val == 10 ? "hartleys" : nothing
    print(io, typeof(scores), " log base=", base_val)
    if unit !== nothing
        print(io, " units=", unit)
    end
    print(io, " : ")
    print(io, "\ntable : ")
    show(io, MIME"text/plain"(), gettable(scores))
end

function Base.show(io::IO, ::MIME"text/plain", pfm::PositionFrequencyMatrix)
    print(io, typeof(pfm), " : ")
    print(io, "\ntable : ")
    show(io, MIME"text/plain"(), gettable(pfm))
end

function Base.show(io::IO, ::MIME"text/plain", ppm::PositionSpecificProbabilityMatrix)
    print(io, typeof(ppm), " : ")
    print(io, "\ntable : ")
    show(io, MIME"text/plain"(), gettable(ppm))
end

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

# Frequencies
# ===========

function _position_frequency_matrix(
    msa::AbstractArray{Residue};
    alphabet::ResidueAlphabet = UngappedAlphabet(),
    weights::WeightTypes = NoClustering(),
    pseudocounts::Pseudocount = NoPseudocount(),
)
    nres = length(alphabet)
    ncols = ncolumns(msa)
    frequencies = Matrix{Float64}(undef, nres, ncols)

    column_table = ContingencyTable(Float64, Val{1}, alphabet)

    for j = 1:ncols
        cleanup!(column_table)
        col_view = @inbounds @view msa[:, j]
        frequencies!(column_table, col_view; weights = weights, pseudocounts = pseudocounts)
        @inbounds frequencies[:, j] .= gettablearray(column_table)
    end

    row_dict = getnamedict(alphabet)
    col_names = columnnames(msa)
    col_dict = OrderedDict{String,Int}(col_names[i] => i for i = 1:ncols)
    NamedArray(frequencies, (row_dict, col_dict), ("Res", "Col"))
end

"""
    position_frequency_matrix(msa::AbstractArray{Residue}; kwargs...)

Compute a position frequency matrix (PFM) from a protein MSA. Each column summarizes how
often each residue appears at that alignment position. Therefore, for column i and residue
a, ``fᵢ(a)`` is the (possibly weighted and pseudocounted) count of that residue in that
column. The `weights` reduce redundancy in highly similar MSAs (so closely related
sequences do not dominate), and `pseudocounts` smooth sparse columns to avoid zero
frequencies. The chosen `alphabet` defines which symbols are counted; for example, an
ungapped alphabet ignores gaps while a gapped alphabet treats the gap as a symbol.

# Keyword Arguments

  - `alphabet = UngappedAlphabet()`: residue alphabet defining the row order.
  - `weights = NoClustering()`: sequence weights to reduce redundancy.
  - `pseudocounts = NoPseudocount()`: smoothing applied before normalization.

# Returns

[`PositionFrequencyMatrix`](@ref) with $(_DOC_PROFILE_SHAPE).
"""
function position_frequency_matrix(
    msa::AbstractArray{Residue};
    alphabet::ResidueAlphabet = UngappedAlphabet(),
    weights::WeightTypes = NoClustering(),
    pseudocounts::Pseudocount = NoPseudocount(),
)
    table = _position_frequency_matrix(
        msa;
        alphabet = alphabet,
        weights = weights,
        pseudocounts = pseudocounts,
    )
    PositionFrequencyMatrix(table, alphabet)
end

# Probabilities
# ============

function _normalize_position_probabilities!(matrix::AbstractMatrix)
    @inbounds for j = 1:size(matrix, 2)
        col = @view matrix[:, j]
        total = sum(col)
        if total == 0.0
            fill!(col, NaN)
        elseif total != 1.0
            col ./= total
        end
    end
    matrix
end

"""
    position_specific_probability_matrix(pfm::PositionFrequencyMatrix)
    position_specific_probability_matrix(msa::AbstractArray{Residue}; kwargs...)

Build a position-specific probability matrix (PSPM). Internally, it converts a PFM
(position frequency matrix) into a PSPM by column-wise normalization using
$(_DOC_PSPM_FORMULA) When built from an MSA, the PFM is computed first using the provided
keyword arguments.

If a column has zero total weight (e.g., all gaps under an ungapped alphabet), that
column is filled with `NaN`.

# Keyword Arguments

$(_DOC_MSA_KWARGS)

# Returns

[`PositionSpecificProbabilityMatrix`](@ref) with $(_DOC_PROFILE_SHAPE).
"""
function position_specific_probability_matrix(pfm::PositionFrequencyMatrix)
    table = deepcopy(pfm.table)
    _normalize_position_probabilities!(getarray(table))
    PositionSpecificProbabilityMatrix(table, pfm.alphabet)
end

function position_specific_probability_matrix(
    msa::AbstractArray{Residue};
    alphabet::ResidueAlphabet = UngappedAlphabet(),
    weights::WeightTypes = NoClustering(),
    pseudocounts::Pseudocount = NoPseudocount(),
)
    table = _position_frequency_matrix(
        msa;
        alphabet = alphabet,
        weights = weights,
        pseudocounts = pseudocounts,
    )
    _normalize_position_probabilities!(getarray(table))
    PositionSpecificProbabilityMatrix(table, alphabet)
end

# Scoring
# =======

"""
    position_specific_scoring_matrix(msa::AbstractArray{Residue}; kwargs...)
    position_specific_scoring_matrix(pfm::PositionFrequencyMatrix; kwargs...)
    position_specific_scoring_matrix(ppm::PositionSpecificProbabilityMatrix; kwargs...)

Compute a position-specific scoring matrix (PSSM) with log-odds scores
$(_DOC_PSSM_FORMULA) The PSSM can be built from an MSA, a position frequency matrix (PFM),
or a position-specific probability matrix (PSPM).

Positive scores indicate enrichment over background; negative scores indicate depletion.
PSSMs are used to score sequences against a profile.

# Keyword Arguments

$(_DOC_MSA_KWARGS)

  - `background = BLOSUM62_Pi`: background distribution `q(a)`; accepts `AbstractArray`,
    `Probabilities` or `ContingencyTable` objects. It is normalized if needed.
  - `base::Number = ℯ`: logarithm base.

Use the keyword argument `base` to change the base of the log. $_DOC_LOG_BASE

Gaps are handled by the chosen `alphabet`: `UngappedAlphabet()` ignores gaps, while
`GappedAlphabet()` includes them. The `background` distribution must match the alphabet
length and is normalized if needed.

Scores are computed directly as log-odds, so IEEE floating-point rules apply for zeros:
pᵢ(a) = 0 ⇒ `-Inf`, q(a) = 0 ⇒ `Inf`, invalid ratios ⇒ `NaN`. Columns with no valid
observations (e.g. all gaps using `UngappedAlphabet()`) are filled with `NaN`.

# Notes

$(_LOG_BASE_UNITS)

# Results

[`PositionSpecificScoringMatrix`](@ref) with $(_DOC_PROFILE_SHAPE).

# Examples

```julia
using MIToS.Information, MIToS.MSA

msa = permutedims(hcat(res"AC", res"AD", res"AE")) # three sequences, two positions
result = position_specific_scoring_matrix(msa)
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
    ppm = position_specific_probability_matrix(
        msa;
        alphabet = alphabet,
        weights = weights,
        pseudocounts = pseudocounts,
    )
    position_specific_scoring_matrix(ppm; background = background, base = base)
end

function position_specific_scoring_matrix(
    ppm::PositionSpecificProbabilityMatrix;
    background = BLOSUM62_Pi,
    base::Number = ℯ,
)
    _check_log_base(base)

    q = _collect_background(background, ppm.alphabet)
    table = deepcopy(ppm.table)
    X = getarray(table)
    nrows = size(X, 1)
    ncols = size(X, 2)
    invlogbase = base === ℯ ? 1.0 : inv(log(base))

    for j = 1:ncols
        for i = 1:nrows
            p = @inbounds X[i, j]
            qi = @inbounds q[i]
            ratio = p / qi
            @inbounds X[i, j] = invlogbase == 1.0 ? log(ratio) : log(ratio) * invlogbase
        end
    end

    PositionSpecificScoringMatrix(table, ppm.alphabet, _convert_base(base))
end

function position_specific_scoring_matrix(
    pfm::PositionFrequencyMatrix;
    background = BLOSUM62_Pi,
    base::Number = ℯ,
)
    ppm = position_specific_probability_matrix(pfm)
    position_specific_scoring_matrix(ppm; background = background, base = base)
end

# Sequence Scoring
# ================

@doc """
    score_sequence(pssm::PositionSpecificScoringMatrix, seq::AbstractVector{Residue})
    score_sequence(pssm::PositionSpecificScoringMatrix, seqs::AbstractMatrix{Residue})
    score_sequence(ppm::PositionSpecificProbabilityMatrix, seq::AbstractVector{Residue}; base = ℯ)
    score_sequence(ppm::PositionSpecificProbabilityMatrix, seqs::AbstractMatrix{Residue}; base = ℯ)

Score sequences against a profile.

The result is a [`ProfileScore`](@ref) with `score`, `kind`, `base`, and `used_positions`
(count of non-gap residues of the sequence that are present in the alphabet). For PSSMs
([`PositionSpecificScoringMatrix`](@ref)), `score` is the sum of per-position log-odds
scores (`kind = :log_odds`). For PSPMs ([`PositionSpecificProbabilityMatrix`](@ref)), 
`score` is the log-likelihood (`kind = :log_likelihood`), computed as the sum of the log
probabilities of the observed residues under the profile. The `base` is the log base used
to compute the scores.

$(_LOG_BASE_UNITS)

Useful transformations of the log-likelihood (`score` when `kind == :log_likelihood`) for 
the sequence under the PSPM (here, `used_positions` is how many positions were scored, and 
`base` is the log base used to compute the `score`):

- The **Average log-likelihood (per position)** is a length-normalized score you can 
  compare across sequences of different lengths: `score / used_positions`
- The **Geometric mean probability (per position)** is a "typical per-position probability" 
  derived from the average log-likelihood: `base^(score / used_positions)`
- The **Likelihood** is the probability of the whole sequence under the PSPM (this may 
  underflow for long sequences, so prefer using `score`): `base^score`

Residues not present in the profile alphabet are skipped. Non-gap residues that are
missing from the alphabet emit one warning per residue/alphabet pair.
""" score_sequence

"""
    ProfileScore(score, kind, base, used_positions)

Container for `score_sequence` results. `score` is a log-odds or log-likelihood sum,
`kind` is `:log_odds` or `:log_likelihood`, `base` is the logarithm base used for the
score, and `used_positions` counts the residues that were present in the profile
alphabet.

$(_LOG_BASE_UNITS)
"""
struct ProfileScore
    score::Float64
    kind::Symbol
    base::Union{Irrational{:ℯ},Float64,Int}
    used_positions::Int
end

function _warn_unknown_residue(res::Residue, alphabet::ResidueAlphabet)
    @warn "Residue $res not in alphabet $(typeof(alphabet)); skipping." maxlog = 1 _id =
        (:score_sequence, res, join(names(alphabet), ','))
end

function _check_sequence_length(psm::AbstractPositionSpecificMatrix, seq)
    ncols = size(psm, 2)
    nres = length(seq)
    if nres != ncols
        throw(ArgumentError("Sequence length $nres doesn't match profile length $ncols."))
    end
    ncols
end

function _sequence_vec(seq::AbstractMatrix{Residue})
    if size(seq, 1) != 1 && size(seq, 2) != 1
        throw(
            ArgumentError(string(
                "A sequence must be represented as a matrix with a singleton dimension, ",
                "but got size $(size(seq)).",)
            ),
        )
    end
    vec(seq)
end

function _alphabet_index(res::Residue, alphabet::ResidueAlphabet)
    if res in alphabet
        return alphabet[res]
    end
    if res != GAP
        _warn_unknown_residue(res, alphabet)
    end
    0
end

function _score_sequence_kernel(
    table::AbstractMatrix{Float64},
    alphabet::ResidueAlphabet,
    seq::AbstractVector{Residue},
    ::Val{log_scores},
    invlogbase::Float64,
) where {log_scores}
    score = 0.0
    used_positions = 0
    @inbounds for i = 1:length(seq)
        idx = _alphabet_index(seq[i], alphabet)
        if idx != 0
            used_positions += 1
            value = table[idx, i]
            if log_scores
                score += log(value) * invlogbase
            else
                score += value
            end
        end
    end
    score, used_positions
end

function _score_sequence(pssm::PositionSpecificScoringMatrix, seq::AbstractVector{Residue})
    _check_sequence_length(pssm, seq)
    alphabet = pssm.alphabet
    table = gettablearray(pssm)
    score, used_positions = _score_sequence_kernel(table, alphabet, seq, Val(false), 1.0)
    ProfileScore(score, :log_odds, pssm.base, used_positions)
end

function _score_sequence(pssm::PositionSpecificScoringMatrix, seqs::AbstractMatrix{Residue})
    _score_sequence(pssm, _sequence_vec(seqs))
end

function _score_sequence(
    ppm::PositionSpecificProbabilityMatrix,
    seq::AbstractVector{Residue};
    base::Number = ℯ,
)
    _check_log_base(base)
    base_val = _convert_base(base)
    _check_sequence_length(ppm, seq)
    alphabet = ppm.alphabet
    table = gettablearray(ppm)
    invlogbase = base_val === ℯ ? 1.0 : inv(log(base_val))
    score, used_positions = _score_sequence_kernel(table, alphabet, seq, Val(true), invlogbase)
    ProfileScore(score, :log_likelihood, base_val, used_positions)
end

function _score_sequence(
    ppm::PositionSpecificProbabilityMatrix,
    seqs::AbstractMatrix{Residue};
    base::Number = ℯ,
)
    _score_sequence(ppm, _sequence_vec(seqs); base = base)
end

score_sequence(pssm::PositionSpecificScoringMatrix, seq::AbstractVector{Residue}) =
    _score_sequence(pssm, seq)

score_sequence(pssm::PositionSpecificScoringMatrix, seqs::AbstractMatrix{Residue}) =
    _score_sequence(pssm, seqs)

score_sequence(
    ppm::PositionSpecificProbabilityMatrix,
    seq::AbstractVector{Residue};
    base::Number = ℯ,
) = _score_sequence(ppm, seq; base = base)

score_sequence(
    ppm::PositionSpecificProbabilityMatrix,
    seqs::AbstractMatrix{Residue};
    base::Number = ℯ,
) = _score_sequence(ppm, seqs; base = base)
