isdefined(Base, :__precompile__) && __precompile__()

module MIToS

  export Utils, MSA, Clustering, Information, PDB, SIFTS, Pfam

  include(joinpath("Utils", "Utils.jl"))
  include(joinpath("MSA", "MSA.jl"))
  include(joinpath("Clustering", "Clustering.jl"))
  include(joinpath("Information", "Information.jl"))
  include(joinpath("PDB", "PDB.jl"))
  include(joinpath("SIFTS", "SIFTS.jl"))
  include(joinpath("Pfam", "Pfam.jl"))

end # module
