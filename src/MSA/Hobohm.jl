# Hobohm I
# ========

"""
Fill `cluster` and `clustersize` matrices. These matrices are assumed to be empty
(only zeroes) and their length is assumed to be equal to the number of sequences
in the alignment (`aln`). `within_cluster` is a predicate function that takes two
sequences and a `threshold` value and returns `true` if they should belong to the
same cluster. `threshold` is passed as the last argument to
`within_cluster`.
"""
function _fill_hobohmI!(
    within_cluster::Function,
    cluster::Vector{Int},
    clustersize::Vector{Int},
    aln::Vector{Vector{Residue}},
    threshold,
)
    cluster_id = 0
    nseq = length(aln)
    @inbounds for i = 1:(nseq-1)
        if cluster[i] == 0
            cluster_id += 1
            cluster[i] = cluster_id
            clustersize[cluster_id] += 1
            ref_seq = aln[i]
            for j = (i+1):nseq
                if cluster[j] == 0 && within_cluster(ref_seq, aln[j], threshold)
                    cluster[j] = cluster_id
                    clustersize[cluster_id] += 1
                end
            end
        end
    end
    @inbounds if cluster[nseq] == 0
        cluster_id += 1
        cluster[nseq] = cluster_id
        clustersize[cluster_id] += 1
    end
    resize!(clustersize, cluster_id)
end

"""
Calculates the weight of each sequence in a cluster. The weight is equal to one divided
by the number of sequences in the cluster.
"""
function _get_sequence_weight(clustersize, cluster)
    nseq = length(cluster)
    sequence_weight = Array{Float64}(undef, nseq)
    for i = 1:nseq
        @inbounds sequence_weight[i] = 1.0 / clustersize[cluster[i]]
    end
    Weights(sequence_weight, Float64(length(clustersize)))
end

"""
`hobohmI(within_cluster, msa, threshold)`

Sequence clustering using the Hobohm I method from Hobohm et al.
`within_cluster` is a predicate function that receives two sequences and
`threshold` and returns `true` when the sequences should be clustered together.
The default behavior is obtained with `percentidentity`.

# References

  - [Hobohm, Uwe, et al. "Selection of representative protein data sets."
    Protein Science 1.3 (1992): 409-417.](@cite 10.1002/pro.5560010313)
"""
function hobohmI(within_cluster::Function, msa::AbstractMatrix{Residue}, threshold)
    aln = getresiduesequences(msa)
    nseq = length(aln)
    cluster = zeros(Int, nseq)
    clustersize = zeros(Int, nseq)
    _fill_hobohmI!(within_cluster, cluster, clustersize, aln, threshold)
    Clusters(clustersize, cluster, _get_sequence_weight(clustersize, cluster))
end

"""
`hobohmI(msa, threshold)`

Convenience method for Hobohm I clustering using `percentidentity` as the
predicate.
"""
hobohmI(msa::AbstractMatrix{Residue}, threshold) = hobohmI(percentidentity, msa, threshold)
