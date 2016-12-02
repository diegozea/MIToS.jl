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

# using Base.Threads

using DataStructures        # OrderedDicts for Annotations
using AutoHashEquals        # Annotations
using NamedArrays           # Col and Seq names, basic sequence/MSA object
using FastaIO               # FastaReader (fast)
using Base.Random           # shuffle!
using PairwiseListMatrices  # Percent Identity Matrices
using Clustering            # Used for sequence clustering: ClusteringResult
using StatsBase             # WeightVec for clustering
using MIToS.Utils

using RecipesBase           # Plots for MSAs

import Clustering: ClusteringResult, nclusters, counts, assignments

export  # Residue
        Residue,
        GAP, XAA,
        @res_str,
        # Alphabet
        ResidueAlphabet,
        GappedAlphabet, UngappedAlphabet, ReducedAlphabet,
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
        # PLM
        sequencepairsmatrix, columnpairsmatrix,
        # Identity
        percentidentity, meanpercentidentity, percentsimilarity,
        # Clusters
        ClusteringResult, # from Clustering.jl
        nclusters, counts, assignments, # from Clustering.jl
        NoClustering, SequenceClusters,
        getweight, nsequences,
        # Hobohm
        hobohmI


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
include("Shuffle.jl")
include("PLM.jl")
include("Identity.jl")
include("Clusters.jl")
include("Hobohm.jl")
# include("Plots.jl")

end
