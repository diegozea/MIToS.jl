# Clusters
# ========

abstract AbstractClusters

immutable NoClustering <: AbstractClusters end

"""
Data structure to represent sequence clusters.
The sequence data itself is not included.
"""
immutable Clusters
    clustersize::Vector{Int}
    sequencecluster::Vector{Int}
    sequenceweight::Vector{Float64}
end

"Get the number of clusters in the set."
nclusters(c::Clusters) = length(c.clustersize)

"Get the weights of all clusters in the set."
getweight(c::Clusters) = c.sequenceweight

"""```getweight(c, i::Int)```

This function returns the weight of the sequence number `i`.
getweight should be defined for any type used for `count!`/`count` in order to use his weigths.
"""
getweight(c::Clusters, seq::Int) = c.sequenceweight[seq]
@inline getweight(weight::NoClustering, i::Int) = 1
getweight(weights::AbstractVector, i::Int) = weights[i]

nsequences(c::Clusters) = length(c.sequencecluster)
