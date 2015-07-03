module PDB

  export PDBResidueIdentifier, Coordinates, PDBAtom, PDBResidue,
  distance, contact,

  getatoms, getresidues


  include("PDBResidues.jl")
  include("PDBMLParser.jl")

end
