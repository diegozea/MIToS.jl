# Clusters
# ========

import MIToS.MSA: nsequences

immutable Clusters
  clustersize::Vector{Int}
  sequencecluster::Vector{Int}
  sequenceweight::Vector{Float64}
end

# Number of clusters
getnclusters(c::Clusters) = length(c.clustersize)

getweight(c::Clusters, seq::Int) = c.sequenceweight[seq]

nsequences(c::Clusters) = length(c.sequencecluster)
