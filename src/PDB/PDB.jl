module PDB

  using LightXML
  using AutoHashEquals
  using DataStructures
  using Formatting
  using MIToS.Utils

  import Base: ==, hash, length, -, norm, dot, angle, cross, vec, any, print, show
  import MIToS.Utils: findobjects, isobject

  export covalentradius, vanderwaalsradius,

  PDBResidueIdentifier, Coordinates, PDBAtom, PDBResidue,
  distance, contact, findheavy, findatoms, findCB, selectbestoccupancy,
  angle,

  ishydrophobic, isaromatic, iscationic, isanionic,
  ishbonddonor, ishbondacceptor, hydrogenbond,
  vanderwaals, vanderwaalsclash, covalent, disulphide,
  aromaticsulphur, pication, aromatic, ionic, hydrophobic,
  stridehydrogenbond, chimerahydrogenbond,

  getpdbmlatoms, downloadpdb,

  getpdbatoms,

  isobject, findobjects, Is, Not

  include("AtomsData.jl")
  include("PDBResidues.jl")
  include("Interaction.jl")
  include("PDBMLParser.jl")
  include("PDBParser.jl")

end
