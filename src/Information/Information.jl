"""
The `Information` module of MIToS defines types and functions useful to calculate
information measures (e.g. *Mutual Information* (MI) and *Entropy*) over a Multiple
Sequence Alignment (MSA). This module was designed to count `Residue`s
(defined in the `MSA` module) in special contingency tables (as fast as possible) and
to derive probabilities from this counts. Also, includes methods for applying corrections
to that tables, e.g. pseudocounts and pseudo frequencies. Finally, `Information` allows to
use this probabilities and counts to estimate information measures and other
frequency based values.

**Features**

- Estimate multi dimensional frequencies and probabilities tables from sequences, MSAs, etc...
- Correction for small number of observations
- Correction for data redundancy on a MSA
- Estimate information measures
- Calculate corrected mutual information between residues

```julia
using MIToS.Information
```
"""
module Information

using MIToS.Utils
using MIToS.MSA
using Serialization         # GaussDCA
using Base.Cartesian        # nloops for ContingencyTables
using NamedArrays           # ContingencyTables have NamedArrays
using DataStructures        # OrderedDicts for NamedArrays
using StatsBase             # entropy
using LinearAlgebra         # normalize
using PairwiseListMatrices

export  # MIToS.MSA
        GappedAlphabet,
        UngappedAlphabet,
        ReducedAlphabet,
        # Pseudocounts
        Pseudocount,
        NoPseudocount,
        AdditiveSmoothing,
        # ContingencyTables
        ContingencyTable,
        Probabilities,
        Counts,
        getcontingencytable,
        getalphabet,
        gettable, gettablearray,
        getmarginals, getmarginalsarray,
        gettotal,
        # update_marginals!,
        apply_pseudocount!,
        delete_dimensions!, delete_dimensions,
        # BLOSUM62
        BLOSUM62_Pi, BLOSUM62_Pij,
        # Pseudofrequencies
        Pseudofrequencies,
        NoPseudofrequencies,
        BLOSUM_Pseudofrequencies,
        apply_pseudofrequencies!,
        # Counters
        count!,
        probabilities,
        probabilities!,
        # Iterations
        mapcolfreq!,
        mapseqfreq!,
        mapcolpairfreq!,
        mapseqpairfreq!,
        cumulative,
        # InformationMeasures
        entropy,
        marginal_entropy,
        kullback_leibler,
        mutual_information,
        normalized_mutual_information,
        # Corrections
        APC!,
        # CorrectedMutualInformation
        buslje09,
        BLMI,
        # Gaps
        gap_union_percentage,
        gap_intersection_percentage,
        pairwisegapfraction,
        # Externals
        gaussdca,
        # Formats from MIToS.MSA
        Raw, Stockholm, FASTA,
        # Imported from Base (and exported for docs)
        normalize,
        normalize!,
        count

include("Pseudocounts.jl")
include("ContingencyTables.jl")
include("BLOSUM62.jl")
include("Pseudofrequencies.jl")
include("Counters.jl")
include("Iterations.jl") # TO DO: Docs
include("InformationMeasures.jl")
include("Corrections.jl")
include("CorrectedMutualInformation.jl")
include("Gaps.jl")
include("Externals.jl")

end
