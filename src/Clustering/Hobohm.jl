# Hobohm I
function hobohmI(aln::Matrix, threshold::Float64)
	(nseq, nres) = Base.size(aln)
	cluster = zeros(Int,nseq)
	weight = zeros(Int,nseq)
	cluster_id = 0
	for i in 1:(nseq-1)
		if cluster[i] == 0
			cluster_id += 1
			cluster[i] = cluster_id
			weight[cluster_id] += 1
			ref_seq = aln[i,:]
			for j in (i+1):nseq
				if cluster[j] == 0 && percentidentity(ref_seq, aln[j,:]) >= threshold
					cluster[j] = cluster_id
					weight[cluster_id] += 1
				end
			end
		end
	end
	if cluster[nseq] == 0
		cluster_id += 1
		cluster[nseq] = cluster_id
		weight[cluster_id] += 1
	end
  resize!(weight, cluster_id)
	Clusters(weight, cluster, Float64[ 1.0 ./ float(weight[cl]) for cl in cluster])
end

hobohmI(msa::MultipleSequenceAlignment, threshold::FloatingPoint) = hobohmI(msa.msa, threshold)
