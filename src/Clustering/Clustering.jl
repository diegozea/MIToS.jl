module Clustering

  using MIToS.MSA

  export Clusters, getnclusters, getweight, nsequences,

  percentidentity, percentidentity2,

  """
  Cluster sequences using Hobohm method 1.
 
  For more information read: 

  1. Protein Sci. 1992 Mar;1(3):409-17.
  Selection of representative protein data sets.
  Hobohm U(1), Scharf M, Schneider R, Sander C.
  """
  hobohmI

  include("Clusters.jl")
  include("Identity.jl")
  include("Hobohm.jl")

end
