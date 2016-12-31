# Counting
# ========

@generated function _temporal_counts!{T,N,A}(counts::ContingencyTable{T,N,A}, weights,
                                  seqs::Vararg{AbstractVector{Residue},N})
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
                       seqs::Vararg{AbstractVector{Residue},N})
    _temporal_counts!(table, weights, seqs...)
    apply_pseudocount!(table, pseudocounts)
    _update!(table)
    table
end

# Default counters
# ================

function _count{N,A <: ResidueAlphabet}(alphabet::A,
                                        weights, pseudocounts,
                                        seqs::Vararg{AbstractVector{Residue},N}
                                        )::ContingencyTable{Float64,N,A}
    table = ContingencyTable(Float64, Val{N}, alphabet)::ContingencyTable{Float64,N,A}
    count!(table, weights, pseudocounts, seqs...)
    table
end

function Base.count{N}(seqs::Vararg{AbstractVector{Residue},N};
                       alphabet::ResidueAlphabet = UngappedAlphabet(),
                       weights = NoClustering(),
                       pseudocounts::Pseudocount = NoPseudocount())
    _count(alphabet, weights, pseudocounts, seqs...)
end

# Probabilities
# =============

function probabilities!{T,N,A}(table::ContingencyTable{T,N,A},
                               weights,
                               pseudocounts::Pseudocount,
                               pseudofrequencies::Pseudofrequencies,
                               seqs::Vararg{AbstractVector{Residue},N})
    count!(table, weights, pseudocounts, seqs...)
    normalize!(table)
    apply_pseudofrequencies!(table, pseudofrequencies)
    table
end

function _probabilities{N,A <: ResidueAlphabet}(alphabet::A,
                                                weights, pseudocounts, pseudofrequencies,
                                                seqs::Vararg{AbstractVector{Residue},N}
                                                )::ContingencyTable{Float64,N,A}
    table = ContingencyTable(Float64, Val{N}, alphabet)::ContingencyTable{Float64,N,A}
    probabilities!(table, weights, pseudocounts, pseudofrequencies, seqs...)
    table
end

function probabilities{N}(seqs::Vararg{AbstractVector{Residue},N};
                          alphabet::ResidueAlphabet = UngappedAlphabet(),
                          weights = NoClustering(),
                          pseudocounts::Pseudocount = NoPseudocount(),
                          pseudofrequencies::Pseudofrequencies = NoPseudofrequencies())
    _probabilities(alphabet, weights, pseudocounts, pseudofrequencies, seqs...)
end
