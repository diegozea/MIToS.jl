module Clustering

  using MIToS.MSA
  using PairwiseListMatrices

  import MIToS.MSA: nsequences

  export AbstractClusters, NoClustering, Clusters, getnclusters, getweight, nsequences,

  percentidentity, meanpercentidentity,

  hobohmI

  include("Clusters.jl")
  include("Identity.jl")
  include("Hobohm.jl")

end
