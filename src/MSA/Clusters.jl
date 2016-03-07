# Clusters
# ========

# Clustering.jl : ClusteringResult

immutable NoClustering <: ClusteringResult end

"""
Data structure to represent sequence clusters.
The sequence data itself is not included.
"""
immutable SequenceClusters <: ClusteringResult
    clustersize::Vector{Int}
    sequencecluster::Vector{Int}
    sequenceweight::Vector{Float64}
end

nsequences(cl::SequenceClusters) = length(cl.sequencecluster)

# Clustering.jl : nclusters, counts, assignments

"Get the number of clusters in `SequenceClusters`."
nclusters(cl::SequenceClusters) = length(cl.clustersize)

"Get sample counts of clusters as a `Vector`. Each `k` value is the number of samples assigned to the k-th cluster."
counts(cl::SequenceClusters) = cl.clustersize

"Get a vector of assignments, where the `i` value is the index/number of the cluster to which the i-th sequence is assigned."
assignments(cl::SequenceClusters) = cl.sequencecluster

#  convert from ClusteringResult of Clustering.jl to SequenceClusters
# -------------------------------------------------------------------

# needs tests
function convert(::Type{SequenceClusters}, cl::ClusteringResult)
    clustersize = counts(cl)
    sequencecluster = assignments(cl)
    sequenceweight = Float64[ 1.0/clustersize[k] for k in sequencecluster ]
    SequenceClusters(clustersize, sequencecluster, sequenceweight)
end

convert(::Type{SequenceClusters}, cl::SequenceClusters) = cl # no-op

# weights
# -------

"Get the weights of all clusters in the set."
getweight(cl::SequenceClusters) = cl.sequenceweight

getweight(cl::ClusteringResult) = getweight(convert(SequenceClusters, cl))

"""```getweight(c, i::Int)```

This function returns the weight of the sequence number `i`.
getweight should be defined for any type used for `count!`/`count` in order to use his weigths.
"""
getweight(c::SequenceClusters, seq::Int) = c.sequenceweight[seq]

@inline getweight(weight::NoClustering, i::Int) = 1

getweight(weights::AbstractVector, i::Int) = weights[i]

