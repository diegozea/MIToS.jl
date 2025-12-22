# Position-Specific Scoring Matrix (PSSM)
# =======================================

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
  - PSPM ([`PositionSpecificProbabilityMatrix`](@ref)): ``pᵢ(a) = fᵢ(a) / ∑ₐ fᵢ(a)``,
    where i is the alignment column, a is a residue symbol in the alphabet, and ``∑ₐ``
    sums over all residues in that alphabet.
  - PSSM ([`PositionSpecificScoreMatrix`](@ref)): ``Sᵢ(a) = log_b(pᵢ(a) / q(a))``, where
    q(a) is the background probability for residue a and b is the log base (units are
    described in the PSSM docstring).

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
stores the frequency matrix in the `table` field. `table` is a `NamedArray` whose rows
follow the provided `alphabet` ordering and whose columns correspond to alignment positions.
"""
struct PositionFrequencyMatrix{A} <: AbstractPositionSpecificMatrix{A}
    table::NamedArray{Float64,2,Array{Float64,2},NTuple{2,OrderedDict{String,Int}}}
    alphabet::A
end

"""
    struct PositionSpecificProbabilityMatrix{A}
        table::NamedArray{Float64,2,Array{Float64,2},NTuple{2,OrderedDict{String,Int}}}
        alphabet::A
    end

Position-specific probability matrix (PSPM). Each column is a probability distribution
over the alphabet at a position: ``pᵢ(a) = fᵢ(a) / ∑ₐ fᵢ(a)``, where i is the alignment
column, a is a residue symbol in the alphabet, and ``∑ₐ`` sums over all residues in that
alphabet. Columns are normalized independently; if a column has zero total weight (for
example, an all-gap column under an ungapped alphabet), that column is filled with `NaN`.
"""
struct PositionSpecificProbabilityMatrix{A} <: AbstractPositionSpecificMatrix{A}
    table::NamedArray{Float64,2,Array{Float64,2},NTuple{2,OrderedDict{String,Int}}}
    alphabet::A
end

"""
    struct PositionSpecificScoreMatrix{A}
        table::NamedArray{Float64,2,Array{Float64,2},NTuple{2,OrderedDict{String,Int}}}
        alphabet::A
        base::Union{Irrational{:ℯ}, Float64, Int}
    end

Position-specific scoring matrix (PSSM) with log-odds scores relative to a background
distribution. Scores are defined as ``Sᵢ(a) = log_b(pᵢ(a) / q(a))``, where i is the
alignment column, a is a residue symbol in the alphabet, pᵢ(a) is the PSPM probability,
q(a) is the background probability for residue a, and b is the log base.

Positive scores indicate enrichment over background; negative scores indicate depletion.
The `base` field controls units: b = 2 yields bits, b = ℯ yields nats, and b = 10 yields
hartleys. IEEE floating-point rules apply for zeros (pᵢ(a) = 0 ⇒ `-Inf`, q(a) = 0 ⇒
`Inf`, invalid ratios ⇒ `NaN`).

To score an aligned sequence x of length n_positions, you can sum per-position scores:
``score(x) = ∑ᵢ Sᵢ(xᵢ)``, where xᵢ is the residue observed at position i.

`table` is a `NamedArray` whose rows follow the provided `alphabet` ordering and whose
columns match the alignment positions.
"""
struct PositionSpecificScoreMatrix{A} <: AbstractPositionSpecificMatrix{A}
    table::NamedArray{Float64,2,Array{Float64,2},NTuple{2,OrderedDict{String,Int}}}
    alphabet::A
    base::Union{Irrational{:ℯ},Float64,Int}
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

function Base.show(io::IO, ::MIME"text/plain", scores::PositionSpecificScoreMatrix)
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

[`PositionFrequencyMatrix`](@ref) with rows in alphabet order and columns matching
alignment positions.
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
``pᵢ(a) = fᵢ(a) / ∑ₐ fᵢ(a)``, where i is the alignment column, a is a residue symbol in
the alphabet, and ``∑ₐ`` sums over all residues in that alphabet. When built from an MSA,
the PFM is computed first using the provided keyword arguments.

If a column has zero total weight (e.g., all gaps under an ungapped alphabet), that
column is filled with `NaN`.

# Keyword Arguments

  - `alphabet = UngappedAlphabet()`: residue alphabet defining the row order.
  - `weights = NoClustering()`: sequence weights to reduce redundancy.
  - `pseudocounts = NoPseudocount()`: smoothing applied before normalization.

# Returns

[`PositionSpecificProbabilityMatrix`](@ref) with rows in alphabet order and columns
matching alignment positions.
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
``Sᵢ(a) = log_b(pᵢ(a) / q(a))``. Here i is the alignment column, a is a residue symbol in
the alphabet, pᵢ(a) is the PSPM probability at that column, q(a) is the background
probability for residue a, and b is the log base. The PSSM can be built from an MSA, a
position frequency matrix (PFM), or a position-specific probability matrix (PSPM).

Positive scores indicate enrichment over background; negative scores indicate depletion.
PSSMs are used to score sequences against a profile.

# Keyword Arguments

  - `alphabet = UngappedAlphabet()`: residue alphabet; score rows follow this order.
  - `weights = NoClustering()`: sequence weights used during probability estimation.
  - `pseudocounts = NoPseudocount()`: pseudocounts applied before normalization.
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

The base controls units: b = 2 yields bits, b = ℯ yields nats, and b = 10 yields
hartleys.

# Results

[`PositionSpecificScoreMatrix`](@ref) with rows in alphabet order and columns matching
alignment positions.

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
    if !(base > 0) || base == 1 # this also catches NaN
        throw(
            ArgumentError(
                "The logarithm base must be positive and different from 1 (base=$base).",
            ),
        )
    end

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

    PositionSpecificScoreMatrix(table, ppm.alphabet, _convert_base(base))
end

function position_specific_scoring_matrix(
    pfm::PositionFrequencyMatrix;
    background = BLOSUM62_Pi,
    base::Number = ℯ,
)
    ppm = position_specific_probability_matrix(pfm)
    position_specific_scoring_matrix(ppm; background = background, base = base)
end
