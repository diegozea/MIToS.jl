module MIToS

  include(joinpath("MSA", "MSA.jl"))
  include(joinpath("Clustering", "Clustering.jl"))
  include(joinpath("PDB", "PDB.jl"))
  include(joinpath("SIFTS", "SIFTS.jl"))

end # module
