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

"""```getweight(c, i::Int)```

This function returns the weight of the sequence number `i`.
getweight should be defined for any type used for `count!`/`count` in order to use his weigths.
"""
getweight(c::Clusters, seq::Int) = c.sequenceweight[seq]
@inline getweight(weight::Real, i::Int) = weight
getweight(weights::AbstractVector, i::Int) = c[i]

nsequences(c::Clusters) = length(c.sequencecluster)
