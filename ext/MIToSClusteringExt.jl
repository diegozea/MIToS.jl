module MIToSClusteringExt

using MIToS.MSA
using Clustering

using StatsBase

"""
Get the number of clusters in a `Clusters` object.
"""
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

# Convert from ClusteringResult to Clusters

# needs tests
function Base.convert(::Type{Clusters}, cl::ClusteringResult)
    clustersize = counts(cl)
    clusters = assignments(cl)
    weights = StatsBase.Weights(
        Float64[1.0 / clustersize[k] for k in clusters],
        Float64(length(clustersize)),
    )
    Clusters(clustersize, clusters, weights)
end

# getweight(cl::ClusteringResult) = getweight(convert(Clusters, cl))
# getweight(cl::ClusteringResult, seq::Int) = getweight(convert(Clusters, cl), seq)

end
