# Clusters
# ========

abstract type AbstractCluster end

"""
Use `NoClustering()` to avoid the use of clustering where a `Clusters` type is needed.
"""
struct NoClustering <: AbstractCluster end

"""
Data structure to represent sequence clusters. The sequence data itself is not included.
"""
@auto_hash_equals struct Clusters <: AbstractCluster
    clustersize::Vector{Int}
    clusters::Vector{Int}
    weights::StatsBase.Weights{Float64,Float64,Array{Float64,1}}
end

nelements(cl::Clusters) = length(cl.clusters)

@inline Base.convert(::Type{Clusters}, cl::Clusters) = cl # no-op

# weights
# -------

"""
The `WeightTypes` type is the same as `Union{Weights,NoClustering,Clusters}`. This type is
used to represent weights. Most of the functions taking the `weights` kerword argument in
the `Information` module accept instances of `WeightTypes`.
"""
const WeightTypes = Union{Weights,NoClustering,Clusters}

"""
`getweight(c[, i::Int])`

This function returns the weight of the sequence number `i`. getweight should be defined for
any type used for `frequencies!`/`frequencies` in order to use his weigths. If `i` isn't
used, this function returns a vector with the weight of each sequence.
"""
@inline getweight(weight::NoClustering, seq::Int) = 1.0

getweight(cl::Clusters) = cl.weights

getweight(cl::Clusters, seq::Int) = cl.weights[seq]

@inline getweight(cl::Weights, i::Int) = cl[i]


