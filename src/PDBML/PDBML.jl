module PDBML

  export PDBResidueIdentifier, Coordinates, PDBAtom, PDBResidue,
  distance, contact, getoccupancy,

  getatoms, getresidues


  include("PDBResidues.jl")
  include("PDBMLParser.jl")

end
