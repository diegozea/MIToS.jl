# Hobohm I
# ========

"""
Fill `cluster` and `clustersize` matrices.
These matrices are assumed to be empty (only zeroes) and their length is assumed to be equal to the number
of sequences in the alignment (`aln`).
`threshold` is the minimum identity value between two sequences to be in the same cluster.
"""
function _fill_hobohmI!(cluster::Vector{Int}, clustersize::Vector{Int}, aln::Vector{Vector{Residue}}, threshold)
    cluster_id = 0
    nseq = length(aln)
    for i in 1:(nseq-1)
        if cluster[i] == 0
            cluster_id += 1
            cluster[i] = cluster_id
            clustersize[cluster_id] += 1
            ref_seq = aln[i]
            for j in (i+1):nseq
                if cluster[j] == 0 && percentidentity(ref_seq, aln[j], threshold)
                    cluster[j] = cluster_id
                    clustersize[cluster_id] += 1
                end
            end
        end
    end
    if cluster[nseq] == 0
        cluster_id += 1
        cluster[nseq] = cluster_id
        clustersize[cluster_id] += 1
    end
    resize!(clustersize, cluster_id)
end

"""
Calculates the weight of each sequence in a cluster.
The weight is equal to one divided by the number of sequences in the cluster.
"""
function _get_sequence_weight(clustersize, cluster)
    nseq = length(cluster)
    sequence_weight = Array(Float64, nseq)
    for i in 1:nseq
        @inbounds sequence_weight[i] = 1.0 / clustersize[cluster[i]]
    end
    sequence_weight
end

"""
Sequence clustering using the Hobohm I method from Hobohm et. al. 1992.

*Hobohm, Uwe, et al. "Selection of representative protein data sets." Protein Science 1.3 (1992): 409-417.*
"""
function hobohmI(msa::AbstractMatrix{Residue}, threshold)
    aln = getresiduesequences(msa)
    nseq = length(aln)
    cluster = zeros(Int,nseq)
    clustersize = zeros(Int,nseq)
    _fill_hobohmI!(cluster, clustersize, aln, threshold)
    SequenceClusters(clustersize, cluster, _get_sequence_weight(clustersize, cluster))
end
