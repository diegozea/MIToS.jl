"""
The MSA module of MIToS has utilities for working with Multiple Sequence Alignments of protein Sequences (MSA).

**Features**

  - Read and write MSAs in `Stockholm`, `FASTA`, `A3M`, `PIR`, or `Raw` format
  - Handle MSA annotations
  - Edit the MSA, e.g. delete columns or sequences, change sequence order, shuffling...
  - Keep track of positions and annotations after modifications on the MSA
  - Describe a MSA, e.g. mean percent identity, sequence coverage, gap percentage...

```julia
using MIToS.MSA
```
"""
module MSA

using OrderedCollections        # OrderedDicts for Annotations
using AutoHashEquals        # Annotations, Clusters
using NamedArrays           # Col and Seq names, basic sequence/MSA object
using FastaIO               # FastaReader (fast)
using Random                # default_rng, shuffle!, rand, Sampler, randstring
using Dates                 # Dates.now()
using PairwiseListMatrices  # Percent Identity Matrices
using Clustering            # Used for sequence clustering: ClusteringResult
using StatsBase             # Weights for clustering
using RecipesBase           # Plots for MSAs
using TranscodingStreams    # To solve MethodError seek(::TranscodingStream, ::Int)
using MIToS.Utils

import Clustering: ClusteringResult, nclusters, counts, assignments
import Markdown: @md_str # for docstrings

export  # Residue
    Residue,
    GAP,
    XAA,
    @res_str,
    # Alphabet
    ResidueAlphabet,
    GappedAlphabet,
    UngappedAlphabet,
    ReducedAlphabet,
    getnamedict,
    # ThreeLetters
    residue2three,
    three2residue,
    # Annotations
    Annotations,
    # filtersequences!,
    ncolumns,
    filtercolumns!,
    getannotfile,
    getannotcolumn,
    getannotsequence,
    getannotresidue,
    setannotfile!,
    setannotcolumn!,
    setannotsequence!,
    setannotresidue!,
    annotate_modification!,
    delete_annotated_modifications!,
    printmodifications,
    # MultipleSequenceAlignment
    AbstractAlignedObject,
    AbstractMultipleSequenceAlignment,
    AbstractAlignedSequence,
    MultipleSequenceAlignment,
    AnnotatedMultipleSequenceAlignment,
    AlignedSequence,
    AnnotatedAlignedSequence,
    AnnotatedAlignedObject,
    UnannotatedAlignedObject,
    NamedResidueMatrix,
    annotations,
    namedmatrix,
    nsequences,
    getresidues,
    getsequence,
    getresiduesequences,
    stringsequence,
    getcolumnmapping,
    getsequencemapping,
    sequencenames,
    columnnames,
    sequencename_iterator,
    columnname_iterator,
    rename_sequences!,
    rename_sequences, # TODO: rename_columns!
    # MSAStats
    gapfraction,
    residuefraction,
    coverage,
    columngapfraction,
    # MSAEditing
    filtersequences,
    filtersequences!,
    filtercolumns,
    filtercolumns!,
    swapsequences!,
    setreference!,
    adjustreference,
    adjustreference!,
    gapstrip,
    gapstrip!,
    # GetIndex
    sequence_index,
    column_index,
    # GeneralParserMethods
    SequenceFormat,
    MSAFormat,
    deletefullgapcolumns,
    deletefullgapcolumns!,
    # Raw
    Raw,
    # Stockholm
    Stockholm,
    # FASTA
    FASTA,
    # NBRF/PIR
    PIR,
    # A3M
    A3M,
    A2M,
    # Shuffle
    shuffle_msa!,
    shuffle_msa,
    # PLM
    sequencepairsmatrix,
    columnpairsmatrix,
    # Identity
    percentidentity,
    meanpercentidentity,
    percentsimilarity,
    # Clusters
    WeightTypes,
    NoClustering,
    Clusters,
    ClusteringResult, # from Clustering.jl
    nclusters,
    counts,
    assignments, # from Clustering.jl
    getweight,
    nelements,
    # Hobohm
    hobohmI,
    # Imported from Utils
    read_file,
    parse_file,
    write_file,
    print_file,
    # Imported from Base or StdLib (and exported for docs)
    names,
    parse,
    isvalid,
    rand,
    shuffle,
    shuffle!,
    # Concatenation
    gethcatmapping

include("Residues.jl")
include("Alphabet.jl")
include("ThreeLetters.jl")
include("Annotations.jl")
include("MultipleSequenceAlignment.jl")
include("MSAEditing.jl")
include("GetIndex.jl")
include("MSAStats.jl")
include("GeneralParserMethods.jl")
include("Raw.jl")
include("Stockholm.jl")
include("FASTA.jl")
include("PIR.jl")
include("A3M.jl")
include("Shuffle.jl")
include("PLM.jl")
include("Identity.jl")
include("Clusters.jl")
include("Hobohm.jl")
include("Plots.jl")
include("Concatenation.jl")

end
