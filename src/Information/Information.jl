"""
The `Information` module of MIToS defines types and functions useful to calculate
information measures (e.g. *Mutual Information* (MI) and *Entropy*) over a Multiple Sequence Alignment (MSA).
This module was designed to count `Residue`s (defined in the `MSA` module)
in special contingency tables (as fast as possible) and to derive probabilities from this counts.
Also, includes methods for applying corrections to that tables, e.g. pseudocounts and pseudo frequencies.
Finally, `Information` allows to use this probabilities and counts to estimate information measures and other frequency based values.

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

using MIToS: Utils, MSA
using Base.Cartesian        # nloops for ContingencyTables
using NamedArrays           # ContingencyTables have NamedArrays
using PairwiseListMatrices


export  # MIToS.MSA
        GappedAlphabet,
        UngappedAlphabet,
        ReducedAlphabet,
        # Pseudocounts
        Pseudocount,
        AdditiveSmoothing,
        # ContingencyTables
        ContingencyTable,
        get_alphabet,
        get_table,
        get_marginals,
        get_total,
        update_marginals!,
        apply_pseudocount!,
        count!

# BLOSUM62_Pi, BLOSUM62_Pij,
#
# SequenceWeights, Pseudocount, AdditiveSmoothing,
# ResidueContingencyTables, ResidueCount, ResidueProbability,
# nresidues, update!, apply_pseudocount!, count!, normalize!,
# blosum_pseudofrequencies!, apply_pseudofrequencies!, probabilities,
# delete_dimensions!, delete_dimensions,
#
# AbstractMeasure, SymmetricMeasure,
# Entropy, KullbackLeibler,
# MutualInformation, MutualInformationOverEntropy,
# estimate, estimate_on_marginal,
# GapUnionPercentage, GapIntersectionPercentage,
#
# estimateincolumns, estimateinsequences, cumulative,
#
# APC!,
#
# buslje09, BLMI, pairwisegapfraction,
#
# # Formats from MIToS.MSA:
# Raw, Stockholm, FASTA,
#
# # Externals
# gaussdca

include("Pseudocounts.jl")
include("ContingencyTables.jl")

# include("BLOSUM62.jl")
# include("Probabilities.jl")
# include("InformationMeasures.jl")
# include("Iterations.jl")
# include("Corrections.jl")
# include("CorrectedMutualInformation.jl")
# include("Externals.jl")

end
