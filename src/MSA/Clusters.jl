# Clusters
# ========

# Clustering.jl : ClusteringResult

immutable NoClustering <: ClusteringResult end

"Data structure to represent sequence clusters. The sequence data itself is not included."
immutable SequenceClusters <: ClusteringResult
    clustersize::Vector{Int}
    sequencecluster::Vector{Int}
    sequenceweight::StatsBase.WeightVec{Float64,Array{Float64,1}}
end

nsequences(cl::SequenceClusters) = length(cl.sequencecluster)

# Clustering.jl : nclusters, counts, assignments

"Get the number of clusters in a `SequenceClusters` object."
Clustering.nclusters(cl::SequenceClusters) = length(cl.clustersize)

"""
Get sample counts of clusters as a `Vector`. Each `k` value is the number of samples
assigned to the k-th cluster.
"""
Clustering.counts(cl::SequenceClusters) = cl.clustersize

"""
Get a vector of assignments, where the `i` value is the index/number of the cluster to
which the i-th sequence is assigned.
"""
Clustering.assignments(cl::SequenceClusters) = cl.sequencecluster

#  convert from ClusteringResult of Clustering.jl to SequenceClusters
# -------------------------------------------------------------------

# needs tests
function Base.convert(::Type{SequenceClusters}, cl::ClusteringResult)
    clustersize = counts(cl)
    sequencecluster = assignments(cl)
    sequenceweight = WeightVec(Float64[1.0/clustersize[k] for k in sequencecluster],
        Float64(length(clustersize)))
    SequenceClusters(clustersize, sequencecluster, sequenceweight)
end

@inline Base.convert(::Type{SequenceClusters}, cl::SequenceClusters) = cl # no-op

# weights
# -------

"""
`getweight(c[, i::Int])`

This function returns the weight of the sequence number `i`. getweight should be defined for
any type used for `count!`/`count` in order to use his weigths. If `i` isn't used, this
function returns a vector with the weight of each sequence.
"""
@inline getweight(weight::NoClustering, seq::Int) = 1.0

getweight(cl::SequenceClusters) = cl.sequenceweight

getweight(cl::SequenceClusters, seq::Int) = cl.sequenceweight[seq]

getweight(cl::ClusteringResult) = getweight(convert(SequenceClusters, cl))

getweight(cl::ClusteringResult, seq::Int) = getweight(convert(SequenceClusters, cl), seq)
