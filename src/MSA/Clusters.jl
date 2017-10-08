# Clusters
# ========

# Clustering.jl : ClusteringResult

"Use `NoClustering()` to avoid the use of clustering where a `Clusters` type is needed."
struct NoClustering <: ClusteringResult end

"Data structure to represent sequence clusters. The sequence data itself is not included."
@auto_hash_equals struct Clusters <: ClusteringResult
    clustersize::Vector{Int}
    clusters::Vector{Int}
    weights::StatsBase.Weights{Float64, Float64, Array{Float64,1}}
end

nelements(cl::Clusters) = length(cl.clusters)

# Clustering.jl : nclusters, counts, assignments

"Get the number of clusters in a `Clusters` object."
Clustering.nclusters(cl::Clusters) = length(cl.clustersize)

"""
Get sample counts of clusters as a `Vector`. Each `k` value is the number of samples
assigned to the k-th cluster.
"""
Clustering.counts(cl::Clusters) = cl.clustersize

"""
Get a vector of assignments, where the `i` value is the index/number of the cluster to
which the i-th sequence is assigned.
"""
Clustering.assignments(cl::Clusters) = cl.clusters

#  convert from ClusteringResult of Clustering.jl to Clusters
# -------------------------------------------------------------------

# needs tests
function Base.convert(::Type{Clusters}, cl::ClusteringResult)
    clustersize = counts(cl)
    clusters = assignments(cl)
    weights = Weights(Float64[1.0/clustersize[k] for k in clusters],
        Float64(length(clustersize)))
    Clusters(clustersize, clusters, weights)
end

@inline Base.convert(::Type{Clusters}, cl::Clusters) = cl # no-op

# weights
# -------

"""
`getweight(c[, i::Int])`

This function returns the weight of the sequence number `i`. getweight should be defined for
any type used for `count!`/`count` in order to use his weigths. If `i` isn't used, this
function returns a vector with the weight of each sequence.
"""
@inline getweight(weight::NoClustering, seq::Int) = 1.0

getweight(cl::Clusters) = cl.weights

getweight(cl::Clusters, seq::Int) = cl.weights[seq]

@inline getweight(cl::Weights, i::Int) = cl[i]

# getweight(cl::ClusteringResult) = getweight(convert(Clusters, cl))
# getweight(cl::ClusteringResult, seq::Int) = getweight(convert(Clusters, cl), seq)
