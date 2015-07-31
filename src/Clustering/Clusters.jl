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

getweight(c::Clusters) = c.sequenceweight
getweight(c::Clusters, seq::Int) = c.sequenceweight[seq]
@inline getweight(weight::Real, i) = weight

nsequences(c::Clusters) = length(c.sequencecluster)
