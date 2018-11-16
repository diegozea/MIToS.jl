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

using RecipesBase           # Plots for PDB Residues
using AutoHashEquals
using StaticArrays
using DataStructures
using PairwiseListMatrices
using NamedArrays
using LightXML
using LinearAlgebra
using Statistics            # mean
using MIToS.Utils
using Formatting

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
        @residuesdict,
        atoms,
        @atoms,
        findheavy,
        findatoms,
        findCB,
        selectbestoccupancy,
        bestoccupancy,
        residuepairsmatrix,
        proximitymean,
        # AtomsData
        covalentradius,
        vanderwaalsradius,
        check_atoms_for_interactions,
        # Interaction
        ishydrophobic,
        isaromatic,
        iscationic,
        isanionic,
        ishbonddonor,
        ishbondacceptor,
        hydrogenbond,
        vanderwaals,
        vanderwaalsclash,
        covalent,
        disulphide,
        aromaticsulphur,
        pication,
        aromatic,
        ionic,
        hydrophobic,
        stridehydrogenbond,
        chimerahydrogenbond,
        # PDBParser
        PDBFile,
        # PDBMLParser
        PDBML,
        downloadpdb,
        PDBMLHeader,
        downloadpdbheader,
        getpdbdescription,
        # Kabsch
        kabsch,
        center!,
        rmsd,
        getCA,
        CAmatrix,
        coordinatesmatrix,
        change_coordinates,
        centeredcoordinates,
        centeredresidues,
        superimpose,
        mean_coordinates,
        rmsf,
        # MIToS.Utils
        All,
        # Imported from Base (and exported for docs)
        any,
        parse,
        print,
        angle


include("PDBResidues.jl")
include("AtomsData.jl")
include("Interaction.jl")
include("PDBParser.jl")
include("PDBMLParser.jl")
include("Kabsch.jl")
include("Plots.jl")

end
