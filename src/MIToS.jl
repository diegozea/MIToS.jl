module MIToS

export Utils, MSA, Information, PDB, SIFTS, Pfam

include("Utils/Utils.jl")
include("MSA/MSA.jl")
include("Information/Information.jl")
include("PDB/PDB.jl")
include("SIFTS/SIFTS.jl")
include("Pfam/Pfam.jl")

end # module
