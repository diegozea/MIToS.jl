module Clustering

  using MIToS.MSA

  export Clusters, getnclusters, getweight, nsequences,

  percentidentity,

  hobohmI

  include("Clusters.jl")
  include("Identity.jl")
  include("Hobohm.jl")

end
