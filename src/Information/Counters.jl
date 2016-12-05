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

function count!{T,N,A}(table::ContingencyTable{T,N,A}, weights, seqs::AbstractVector{Residue}...)
    _temporal_counts!(table, weights, seqs...)
    _update!(table)
end

# Default counters
# ================

function count!{T,N,A}(table::ContingencyTable{T,N,A}, seqs::AbstractVector{Residue}...)
    count!(table, NoClustering(), seqs...)
end

# Probabilities
# =============

