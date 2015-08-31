isdefined(Base, :__precompile__) && __precompile__(false)

module MIToS

  export Utils, MSA, Clustering, Information, PDB, SIFTS

  include(joinpath("Utils", "Utils.jl"))
  include(joinpath("MSA", "MSA.jl"))
  include(joinpath("Clustering", "Clustering.jl"))
  include(joinpath("Information", "Information.jl"))
  include(joinpath("PDB", "PDB.jl"))
  include(joinpath("SIFTS", "SIFTS.jl"))

end # module
