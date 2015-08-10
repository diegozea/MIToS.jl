# Hobohm I
# ========

function _fill_hobohmI!(cluster::Vector{Int}, clustersize::Vector{Int}, aln::Vector{Vector{Residue}}, threshold::Float64)
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

function _get_sequence_weight(clustersize, cluster)
  nseq = length(cluster)
  sequence_weight = Array(Float64, nseq)
  for i in 1:nseq
    @inbounds sequence_weight[i] = 1.0 / clustersize[cluster[i]]
  end
  sequence_weight
end

function hobohmI(msa::MultipleSequenceAlignment, threshold::Float64)
  aln = getrawsequences(msa)
  nseq = length(aln)
  cluster = zeros(Int,nseq)
  clustersize = zeros(Int,nseq)
  _fill_hobohmI!(cluster, clustersize, aln, threshold)
  Clusters(clustersize, cluster, _get_sequence_weight(clustersize, cluster))
end
