# Counting
# ========

@generated function _temporal_counts!{T,N,A}(counts::ContingencyTable{T,N,A}, weights,
                                  seqs::AbstractVector{Residue}...)
    quote
        @assert N == length(seqs) "Number of residue arrays and table dimension doesn't match."
        # seq_1 = seqs[1]
        # seq_2 = ...
        @nextract $N seq d -> seqs[d]
        len = length(seq_1)
        # @assert len == length(seq_1) "Residue arrays have different lengths"
        # @assert len == length(seq_2) ...
        @nexprs $N d -> @assert len == length(seq_d) "Residue arrays have different lengths."
        if isa(weights, AbstractArray)
            @assert len == length(weights) "Residue array and weights sizes doesn't match."
        end
        _cleanup_temporal!(counts)
        temporal = counts.temporal
        @inbounds @simd for index in 1:length(seq_1)
            # temporal[Int(seq_1[index]), Int(seq_2... += getweight(weights, index)
            @nref($N, temporal, d -> Int(seq_d[index])) += getweight(weights, index)
        end
        counts
    end
end

function count!{T,N,A}(table::ContingencyTable{T,N,A},
                       weights,
                       pseudocounts::Pseudocount,
                       seqs::AbstractVector{Residue}...)
    _temporal_counts!(table, weights, seqs...)
    apply_pseudocount!(table, pseudocounts)
    _update!(table)
    table
end

# Default counters
# ================

function count(seqs::AbstractVector{Residue}...;
               alphabet::ResidueAlphabet = UngappedAlphabet(),
               weights = NoClustering(),
               pseudocounts::Pseudocount = NoPseudocount())
    table = ContingencyTable(Float64, length(seqs), alphabet)
    count!(table, weights, pseudocounts, seqs...)
end

# Probabilities
# =============

function probabilities!{T,N,A}(table::ContingencyTable{T,N,A},
                               weights,
                               pseudocounts::Pseudocount,
                               pseudofrequencies::Pseudofrequencies,
                               seqs::AbstractVector{Residue}...)
    count!(table, weights, pseudocounts, seqs...)
    normalize!(table)
    apply_pseudofrequencies!(table, pseudofrequencies)
    table
end

function probabilities(seqs::AbstractVector{Residue}...;
                       alphabet::ResidueAlphabet = UngappedAlphabet(),
                       weights = NoClustering(),
                       pseudocounts::Pseudocount = NoPseudocount(),
                       pseudofrequencies::Pseudofrequencies = NoPseudofrequencies())
    table = ContingencyTable(Float64, length(seqs), alphabet)
    probabilities!(table, weights, pseudocounts, pseudofrequencies, seqs...)
end
