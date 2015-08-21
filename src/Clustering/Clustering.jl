module Clustering

  using MIToS.MSA

  export Clusters, getnclusters, getweight, nsequences,

  percentidentity, percentidentity2

  hobohmI

  include("Clusters.jl")
  include("Identity.jl")
  include("Hobohm.jl")

end
