module MIToS

  include(joinpath("MSA", "MSA.jl"))
  include(joinpath("Clustering", "Clustering.jl"))
  include(joinpath("PDBML", "PDBML.jl"))
  include(joinpath("SIFTS", "SIFTS.jl"))

end # module
