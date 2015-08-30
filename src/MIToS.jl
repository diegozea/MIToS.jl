VERSION >= v"0.4.0-dev+6521" && __precompile__()

module MIToS

  export Utils, MSA, Clustering, Information, PDB, SIFTS

  include(joinpath("Utils", "Utils.jl"))
  include(joinpath("MSA", "MSA.jl"))
  include(joinpath("Clustering", "Clustering.jl"))
  include(joinpath("Information", "Information.jl"))
  include(joinpath("PDB", "PDB.jl"))
  include(joinpath("SIFTS", "SIFTS.jl"))

end # module
