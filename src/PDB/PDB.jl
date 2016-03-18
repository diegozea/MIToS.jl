"""
The module `PDB` defines types and methods to work with protein structures inside Julia.
It is useful to link structural and sequential information, and
needed for measure the predictive performance at protein contact prediction of mutual information scores.

**Features**

- Read and parse PDF and PDBML files
- Calculate distance and contacts between atoms or residues
- Determine interaction between residues

```julia

using PDB
```
"""
module PDB

  using LightXML
  using AutoHashEquals
  using DataStructures
  using Formatting
  using MIToS.Utils
  # using FixedSizeArrays

  import Base: ==, hash, length, size, -, +, ./, norm, dot, angle, cross, vec, any, print, show, parse
  import MIToS.Utils: findobjects, isobject

  export covalentradius, vanderwaalsradius, check_atoms_for_interactions,

  PDBResidueIdentifier, Coordinates, PDBAtom, PDBResidue,
  distance, contact, findheavy, findatoms, findCB, selectbestoccupancy, bestoccupancy!,
  angle,

  ishydrophobic, isaromatic, iscationic, isanionic,
  ishbonddonor, ishbondacceptor, hydrogenbond,
  vanderwaals, vanderwaalsclash, covalent, disulphide,
  aromaticsulphur, pication, aromatic, ionic, hydrophobic,
  stridehydrogenbond, chimerahydrogenbond,

  PDBFile, PDBML, downloadpdb, getpdbdescription,

  # Mitos.Utils
  isobject, findobjects, Is, Not, In, collectobjects, collectcaptures,

  @residues, residues, @atoms, atoms, @residuesdict, residuesdict

  include("PDBResidues.jl")
  include("AtomsData.jl")
  include("Interaction.jl")
  include("PDBMLParser.jl")
  include("PDBParser.jl")

end
