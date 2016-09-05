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
using IndexedArrays         # IndexedArray for sequence names in MSAs
using PairwiseListMatrices  # Percent Identity Matrices
using Clustering            # Used for sequence clustering: ClusteringResult
using FastaIO               # FastaReader (fast)
using RecipesBase           # Plots for MSAs
using MIToS.Utils

import Base: parse, print, write, convert

import Clustering: ClusteringResult, nclusters, counts, assignments

export Residue, GAP, @res_str, residue2three, three2residue,

Annotations, filtersequences!, filtercolumns!, empty,
getannotfile,  getannotcolumn,  getannotsequence,  getannotresidue,
setannotfile!, setannotcolumn!, setannotsequence!, setannotresidue!,
annotate_modification!, delete_annotated_modifications!, printmodifications,
annotations,

MultipleSequenceAlignment, AnnotatedMultipleSequenceAlignment,
AbstractMultipleSequenceAlignment,
AlignedSequence, AnnotatedAlignedSequence, AbstractAlignedSequence,
getresidues, getsequence, getresiduesequences,
nsequences, ncolumns, gapfraction, residuefraction, coverage,
columngapfraction, setreference!, gapstrip!, adjustreference!, asciisequence,
gapstrip, adjustreference, filtersequences, filtercolumns,
getcolumnmapping, getsequencemapping,

Raw, Stockholm, FASTA,

shuffle_columnwise!, shuffle_sequencewise!, shuffle_residues_sequencewise!,
shuffle_residues_columnwise!,

sequencepairsmatrix, columnpairsmatrix,
columnlabels, sequencelabels,

percentidentity, meanpercentidentity, percentsimilarity,

ClusteringResult, # from Clustering.jl
nclusters, counts, assignments, # from Clustering.jl
NoClustering, SequenceClusters,
getweight, nsequences,

hobohmI,

swap!

include("Residues.jl")
include("Annotations.jl")
include("MultipleSequenceAlignment.jl")
include("GeneralParserMethods.jl")
include("Raw.jl")
include("Stockholm.jl")
include("FASTA.jl")
include("Shuffle.jl")
include("PLM.jl")
include("Identity.jl")
include("Clusters.jl")
include("Hobohm.jl")
include("Plots.jl")

end
