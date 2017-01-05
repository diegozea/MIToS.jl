"""
The module `PDB` defines types and methods to work with protein structures inside Julia.
It is useful to link structural and sequential information, and
needed for measure the predictive performance at protein contact prediction of mutual information scores.

**Features**

- Read and parse PDF and PDBML files
- Calculate distance and contacts between atoms or residues
- Determine interaction between residues

```julia

using MIToS.PDB
```
"""
module PDB

using AutoHashEquals
using FixedSizeArrays
using DataStructures
using Formatting

# using LightXML
# using MIToS.Utils
# using PairwiseListMatrices
# using RecipesBase           # Plots for PDBResidue

# import Base: ==, hash, length, size, -, +, ./, norm, dot, angle, cross, vec, any, print, show, parse
# import MIToS.Utils: findobjects, isobject

export  # PDBResidues
        PDBResidueIdentifier,
        Coordinates,
        PDBAtom,
        PDBResidue,
        squared_distance,
        distance,
        contact,
        isresidue,
        isatom,
        residues,
        @residues,
        residuesdict,
        @residuesdict
        #distance, contact, findheavy, findatoms, findCB, selectbestoccupancy, bestoccupancy,
        #angle, proximitymean

# covalentradius, vanderwaalsradius, check_atoms_for_interactions,
#
#
# ishydrophobic, isaromatic, iscationic, isanionic,
# ishbonddonor, ishbondacceptor, hydrogenbond,
# vanderwaals, vanderwaalsclash, covalent, disulphide,
# aromaticsulphur, pication, aromatic, ionic, hydrophobic,
# stridehydrogenbond, chimerahydrogenbond,
#
# PDBFile, PDBML, downloadpdb, getpdbdescription,
#
# kabsch, center!, rmsd,
# getCA, CAmatrix, coordinatesmatrix, change_coordinates,
# centeredcoordinates, centeredresidues,
# superimpose,
# mean_coordinates, rmsf,
#
# # Mitos.Utils
# isobject, findobjects, Is, Not, In, collectobjects, collectcaptures,
#
# @residues, residues, @atoms, atoms, @residuesdict, residuesdict
#
include("PDBResidues.jl")
# include("AtomsData.jl")
# include("Interaction.jl")
# include("PDBMLParser.jl")
# include("PDBParser.jl")
# include("Kabsch.jl")
# include("Plots.jl")
#
# @deprecate bestoccupancy! bestoccupancy
#
end
