module PDB

  export PDBResidueIdentifier, Coordinates, PDBAtom, PDBResidue,
  distance, contact,

  getpdbmlatoms, getresidues,

  getpdbatoms

  include("PDBResidues.jl")
  include("PDBMLParser.jl")
  include("PDBParser.jl")

end
