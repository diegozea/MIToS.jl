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
using MIToS.Utils

"""
`swap!(ia::IndexedArray, to::Int, from::Int)` interchange/swap the values on the indices `to` and `from` in the `IndexedArray`
"""
function swap!(ia::IndexedArray, to::Int, from::Int)
    previous_id  = ia[to]
    future_id    = ia[from]

    ia.items[to]   = future_id
    ia.items[from] = previous_id

    ia.lookup[previous_id] = from
    ia.lookup[future_id]   = to

    ia
end

import Base: parse, print, write

import Clustering: ClusteringResult, nclusters, counts, assignments

export Residue, GAP, @res_str,
swap!,

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

percentidentity, meanpercentidentity,

ClusteringResult, # from Clustering.jl
nclusters, counts, assignments, # from Clustering.jl
NoClustering, SequenceClusters,
getweight, nsequences,

hobohmI

include("Residues.jl")
include("Annotations.jl")
include("MultipleSequenceAlignment.jl")
include("GeneralParserMethods.jl")
include("Raw.jl")
include("Stockholm.jl")
include("FASTA.jl")
include("Shuffle.jl")
include("Identity.jl")
include("Clusters.jl")
include("Hobohm.jl")

end
