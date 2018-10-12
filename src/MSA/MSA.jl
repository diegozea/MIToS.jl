"""
The MSA module of MIToS has utilities for working with Multiple Sequence Alignments of protein Sequences (MSA).

**Features**

- Read and write MSAs in `Stockholm`, `FASTA` or `Raw` format
- Handle MSA annotations
- Edit the MSA, e.g. delete columns or sequences, change sequence order, shuffling...
- Keep track of positions and annotations after modifications on the MSA
- Describe a MSA, e.g. mean percent identity, sequence coverage, gap percentage...

```julia
using MIToS.MSA
```
"""
module MSA

using DataStructures        # OrderedDicts for Annotations
using AutoHashEquals        # Annotations, Clusters
using NamedArrays           # Col and Seq names, basic sequence/MSA object
using FastaIO               # FastaReader (fast)
using Random                # GLOBAL_RNG, shuffle!, rand, Sampler
using Dates                 # Dates.now()
using PairwiseListMatrices  # Percent Identity Matrices
using Clustering            # Used for sequence clustering: ClusteringResult
using StatsBase             # Weights for clustering
using RecipesBase           # Plots for MSAs
using TranscodingStreams    # To solve MethodError seek(::TranscodingStream, ::Int)
using MIToS.Utils

import Clustering: ClusteringResult, nclusters, counts, assignments

export  # Residue
        Residue,
        GAP, XAA,
        @res_str,
        # Alphabet
        ResidueAlphabet,
        GappedAlphabet, UngappedAlphabet, ReducedAlphabet,
        getnamedict,
        # ThreeLetters
        residue2three, three2residue,
        # Annotations
        Annotations,
        # filtersequences!,
        ncolumns,
        filtercolumns!,
        getannotfile,  getannotcolumn,  getannotsequence,  getannotresidue,
        setannotfile!, setannotcolumn!, setannotsequence!, setannotresidue!,
        annotate_modification!, delete_annotated_modifications!, printmodifications,
        # MultipleSequenceAlignment
        AbstractAlignedObject,
        AbstractMultipleSequenceAlignment,
        AbstractAlignedSequence,
        MultipleSequenceAlignment, AnnotatedMultipleSequenceAlignment,
        AlignedSequence, AnnotatedAlignedSequence,
        AnnotatedAlignedObject, UnannotatedAlignedObject,
        NamedResidueMatrix,
        annotations,
        namedmatrix,
        nsequences,
        getresidues, getsequence, getresiduesequences,
        stringsequence,
        getcolumnmapping, getsequencemapping,
        sequencenames, columnnames, # TO DO: sequencenames!(...)
        # MSAStats
        gapfraction,
        residuefraction,
        coverage,
        columngapfraction,
        # MSAEditing
        filtersequences, filtersequences!,
        filtercolumns, filtercolumns!,
        swapsequences!,
        setreference!,
        adjustreference, adjustreference!,
        gapstrip, gapstrip!,
        # GeneralParserMethods
        deletefullgapcolumns, deletefullgapcolumns!,
        # Raw
        Raw,
        # Stockholm
        Stockholm,
        # FASTA
        FASTA,
        # NBRF/PIR
        PIR,
        # PLM
        sequencepairsmatrix, columnpairsmatrix,
        # Identity
        percentidentity, meanpercentidentity, percentsimilarity,
        # Clusters
        ClusteringResult, # from Clustering.jl
        nclusters, counts, assignments, # from Clustering.jl
        NoClustering, Clusters,
        getweight, nelements,
        # Hobohm
        hobohmI,
        # Imported from Base or StdLib (and exported for docs)
        names,
        parse,
        isvalid,
        rand,
        shuffle,
        shuffle!

include("Residues.jl")
include("Alphabet.jl")
include("ThreeLetters.jl")
include("Annotations.jl")
include("MultipleSequenceAlignment.jl")
include("MSAStats.jl")
include("MSAEditing.jl")
include("GeneralParserMethods.jl")
include("Raw.jl")
include("Stockholm.jl")
include("FASTA.jl")
include("PIR.jl")
include("Shuffle.jl")
include("PLM.jl")
include("Identity.jl")
include("Clusters.jl")
include("Hobohm.jl")
include("Plots.jl")

end
