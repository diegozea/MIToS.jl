module MIToS

  export Utils, MSA, Information, PDB, SIFTS, Pfam

  include(joinpath("Utils", "Utils.jl"))
  include(joinpath("MSA", "MSA.jl"))
  include(joinpath("Information", "Information.jl"))
  include(joinpath("PDB", "PDB.jl"))
  include(joinpath("SIFTS", "SIFTS.jl"))
  include(joinpath("Pfam", "Pfam.jl"))

end # module
