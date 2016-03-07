module Clustering

using MIToS.MSA
using PairwiseListMatrices

import MIToS.MSA: nsequences

export  AbstractClusters, NoClustering, Clusters, nclusters,
getweight, nsequences,

hobohmI

include("Clusters.jl")
include("Hobohm.jl")

end
