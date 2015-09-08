module Clustering

  using MIToS.MSA

  import MIToS.MSA: nsequences

  export AbstractClusters, NoClustering, Clusters, getnclusters, getweight, nsequences,

  percentidentity,

  hobohmI

  include("Clusters.jl")
  include("Identity.jl")
  include("Hobohm.jl")

end
