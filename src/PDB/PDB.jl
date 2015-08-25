module PDB

  using LightXML
  using AutoHashEquals

  import Base: ==, !=, hash, length, -, norm, dot, angle, cross, vec, any

  export covalentradius, vanderwaalsradius,

  PDBResidueIdentifier, Coordinates, PDBAtom, PDBResidue,
  distance, contact, findheavy, findatom, findCB, selectbestoccupancy,
  angle,

  ishydrophobic, isaromatic, iscationic, isanionic,
  ishbonddonor, ishbondacceptor, hydrogenbond,
  vanderwaals, vanderwaalsclash, covalent, disulphide,
  aromaticsulphur, pication, aromatic, ionic, hydrophobic,
  stridehydrogenbond, chimerahydrogenbond,

  getpdbmlatoms, getresidues, downloadpdb,

  getpdbatoms

  include("AtomsData.jl")
  include("PDBResidues.jl")
  include("Interaction.jl")
  include("PDBMLParser.jl")
  include("PDBParser.jl")

end
