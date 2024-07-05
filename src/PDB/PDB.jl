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

import LightXML

using RecipesBase           # Plots for PDB Residues
using AutoHashEquals
using StaticArrays
using OrderedCollections
using PairwiseListMatrices
using NamedArrays
using LinearAlgebra
using Statistics            # mean
using MIToS.Utils
using Format
using JSON3
using Downloads
using Logging

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
    select_residues,
    residues,
    @residues,
    residuesdict,
    @residuesdict,
    select_atoms,
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
    # PDBParser
    PDBFile,
    # PDBMLParser
    PDBML,
    downloadpdb,
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
    read_file,
    parse_file,
    write_file,
    print_file,
    # Sequences
    is_aminoacid,
    modelled_sequences,
    # AlphaFoldDB
    query_alphafolddb,
    download_alphafold_structure,
    # Imported from Base (and exported for docs)
    any,
    parse,
    angle

include("PDBResidues.jl")
include("Sequences.jl")
include("AtomsData.jl")
include("Interaction.jl")
include("PDBParser.jl")
include("PDBMLParser.jl")
include("Kabsch.jl")
include("Plots.jl")
include("AlphaFoldDB.jl")

end
