# Counting
# ========

@generated function _temporal_counts!(
    counts::ContingencyTable{T,N,A},
    weights,
    seqs::Vararg{AbstractVector{Residue},N},
) where {T,N,A}
    quote
        @assert N == length(seqs) "Number of residue arrays and table dimension doesn't match."
        # seq_1 = seqs[1]
        # seq_2 = ...
        @nextract $N seq d -> seqs[d]
        len = length(seq_1)
        # @assert len == length(seq_1) "Residue arrays have different lengths"
        # @assert len == length(seq_2) ...
        @nexprs $N d ->
            @assert len == length(seq_d) "Residue arrays have different lengths."
        if isa(weights, AbstractArray)
            @assert len == length(weights) "Residue array and weights sizes doesn't match."
        end
        _cleanup_temporal!(counts)
        temporal = counts.temporal
        @inbounds @simd for index in eachindex(seq_1)
            # temporal[Int(seq_1[index]), Int(seq_2... += getweight(weights, index)
            @nref($N, temporal, d -> Int(seq_d[index])) += getweight(weights, index)
        end
        counts
    end
end

"""
It populates a `ContingencyTable` (first argument) using the frequencies in the sequences
(last positional arguments). The dimension of the table must match the number of sequences
and all the sequences must have the same length. You must indicate the used weights and
pseudocounts as second and third positional arguments respectively. You can use
`NoPseudofrequencies()` and `NoClustering()` to avoid the use of sequence weighting and
pseudocounts, respectively.

**DEPRECATED**: Use [`frequencies!`](@ref) instead. Note that `frequencies!` defines the weigths and
pseudocounts using keyword arguments instead of positional arguments.
"""
function count!(
    table::ContingencyTable{T,N,A},
    weights,
    pseudocounts::Pseudocount,
    seqs::Vararg{AbstractVector{Residue},N},
) where {T,N,A}
    Base.depwarn(
        "count! using a ContingencyTable or Counts is deprecated. Use frequencies! instead.",
        :count!,
        force = true,
    )
    frequencies!(table, seqs..., weights = weights, pseudocounts = pseudocounts)
end

function count!(table::Counts{T,N,A}, args...) where {T,N,A}
    count!(getcontingencytable(table), args...)
end

# frequencies! is like count! but using keyword arguments
"""
    frequencies!(table, seqs...; weights::WeightTypes, pseudocounts::Pseudocount)

It populates a `ContingencyTable` or `Counts` table (first argument) using the frequencies
in the given sequences (last positional arguments). The dimension of the table must match
the number of sequences and all the sequences must have the same length. You must indicate
the used `weights` and `pseudocounts` as keyword arguments. Those arguments default to
`NoClustering()` and `NoPseudocount()` respectively, to avoid the use of sequence
weighting and pseudocounts.
"""
function frequencies!(
    table::ContingencyTable{T,N,A},
    seqs::Vararg{AbstractVector{Residue},N};
    weights::WeightTypes = NoClustering(),
    pseudocounts::Pseudocount = NoPseudocount(),
) where {T,N,A}
    _temporal_counts!(table, weights, seqs...)
    apply_pseudocount!(table, pseudocounts)
    _update!(table)
    table
end

frequencies!(table::Counts{T,N,A}, args...; kwargs...) where {T,N,A} =
    frequencies!(getcontingencytable(table), args...; kwargs...)

# Default counters
# ================

function _count(
    alphabet::A,
    weights,
    pseudocounts,
    seqs::Vararg{AbstractVector{Residue},N},
)::ContingencyTable{Float64,N,A} where {N,A<:ResidueAlphabet}
    table = ContingencyTable(Float64, Val{N}, alphabet)::ContingencyTable{Float64,N,A}
    count!(table, weights, pseudocounts, seqs...)
    table
end

"""
It returns a `ContingencyTable` wrapped in a `Counts` type with the frequencies of residues
in the sequences that takes as arguments. The dimension of the table is equal to the number
of sequences. You can use the keyword arguments `alphabet`, `weights` and `pseudocounts`
to indicate the alphabet of the table (default to `UngappedAlphabet()`), a clustering
result (default to `NoClustering()`) and the pseudocounts (default to `NoPseudocount()`)
to be used during the estimation of the frequencies.

**DEPRECATED**: Use [`frequencies`](@ref) instead. Note that `frequencies` defines the
alphabet, weigths and pseudocounts using keyword arguments instead of positional arguments.
"""
function Base.count(
    seqs::Vararg{AbstractVector{Residue},N};
    alphabet::ResidueAlphabet = UngappedAlphabet(),
    weights::WeightTypes = NoClustering(),
    pseudocounts::Pseudocount = NoPseudocount(),
) where {N}
    Base.depwarn(
        "`count` on sequences is deprecated in favor of `frequencies`.",
        :count,
        force = true,
    )
    frequencies(
        seqs...;
        alphabet = alphabet,
        weights = weights,
        pseudocounts = pseudocounts,
    )
end


"""
    frequencies(seqs...; alphabet=UngappedAlphabet(), weights=NoClustering(), pseudocounts=NoPseudocount() 

This function returns a `Counts` object wrapping a `ContingencyTable` with the frequencies
of residues in the sequences that takes as arguments. The dimension of the table is equal
to the number of sequences. You can use the keyword arguments `alphabet`, `weights` and
`pseudocounts` to indicate the alphabet of the table, a clustering result and the
pseudocounts to be used during the estimation of the frequencies.
"""
function frequencies(
    seqs::Vararg{AbstractVector{Residue},N};
    alphabet::ResidueAlphabet = UngappedAlphabet(),
    weights::WeightTypes = NoClustering(),
    pseudocounts::Pseudocount = NoPseudocount(),
) where {N}
    Counts(_count(alphabet, weights, pseudocounts, seqs...))
end

# Probabilities
# =============

"""
It populates a `ContingencyTable` (first argument) using the probabilities in the sequences
(last positional arguments). The dimension of the table must match the number of sequences
and all the sequences must have the same length. You must indicate the used weights,
pseudocounts and pseudofrequencies as second, third and fourth positional arguments
respectively. You can use `NoClustering()`, `NoPseudocount()` and `NoPseudofrequencies()`
to avoid the use of sequence weighting, pseudocounts and pseudofrequencies, respectively.
"""
function probabilities!(
    table::ContingencyTable{T,N,A},
    weights,
    pseudocounts::Pseudocount,
    pseudofrequencies::Pseudofrequencies,
    seqs::Vararg{AbstractVector{Residue},N},
) where {T,N,A}
    Base.depwarn(
        "The probabilities! method indicating weights, pseudocounts and pseudofrequencies using positional arguments is deprecated; use keyword arguments instead.",
        :probabilities!,
        force = true,
    )
    probabilities!(
        table,
        seqs...,
        weights = weights,
        pseudocounts = pseudocounts,
        pseudofrequencies = pseudofrequencies,
    )
end

function probabilities!(
    table::ContingencyTable{T,N,A},
    seqs::Vararg{AbstractVector{Residue},N};
    weights::WeightTypes = NoClustering(),
    pseudocounts::Pseudocount = NoPseudocount(),
    pseudofrequencies::Pseudofrequencies = NoPseudofrequencies(),
) where {T,N,A}
    count!(table, weights, pseudocounts, seqs...)
    normalize!(table)
    apply_pseudofrequencies!(table, pseudofrequencies)
    table
end

function probabilities!(table::Probabilities{T,N,A}, args...; kargs...) where {T,N,A}
    probabilities!(getcontingencytable(table), args...; kargs...)
end

# Default probabilities
# =====================

function _probabilities(
    alphabet::A,
    weights,
    pseudocounts,
    pseudofrequencies,
    seqs::Vararg{AbstractVector{Residue},N},
)::ContingencyTable{Float64,N,A} where {N,A<:ResidueAlphabet}
    table = ContingencyTable(Float64, Val{N}, alphabet)::ContingencyTable{Float64,N,A}
    probabilities!(table, weights, pseudocounts, pseudofrequencies, seqs...)
    table
end

"""
It returns a `ContingencyTable` wrapped in a `Probabilities` type with the probabilities of
residues in the sequences that takes as arguments. The dimension of the table is equal to
the number of sequences. You can use the keyword arguments `alphabet`, `weights`,
`pseudocounts` and `pseudofrequencies` to indicate the alphabet of the table
(default to `UngappedAlphabet()`), a clustering result (default to `NoClustering()`),
the pseudocounts (default to `NoPseudocount()`) and the pseudofrequencies
(default to `NoPseudofrequencies()`) to be used during the estimation of the probabilities.
"""
function probabilities(
    seqs::Vararg{AbstractVector{Residue},N};
    alphabet::ResidueAlphabet = UngappedAlphabet(),
    weights::WeightTypes = NoClustering(),
    pseudocounts::Pseudocount = NoPseudocount(),
    pseudofrequencies::Pseudofrequencies = NoPseudofrequencies(),
) where {N}
    Probabilities(
        _probabilities(alphabet, weights, pseudocounts, pseudofrequencies, seqs...),
    )
end
