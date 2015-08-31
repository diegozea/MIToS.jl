module PDB

  using LightXML
  using AutoHashEquals
  using DataStructures
  using Formatting
  using MIToS.Utils

  import Base: ==, hash, length, -, norm, dot, angle, cross, vec, any, print, show, parse
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

  #getpdbmlatoms,
  PDBFile, PDBML, downloadpdb,

  #getpdbatoms,

  # Mitos.Utils
  isobject, findobjects, Is, Not, In

  include("AtomsData.jl")
  include("PDBResidues.jl")
  include("Interaction.jl")
  include("PDBMLParser.jl")
  include("PDBParser.jl")

end
