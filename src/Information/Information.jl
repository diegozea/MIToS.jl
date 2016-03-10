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

  using MIToS.Utils
  using MIToS.MSA
  using PairwiseListMatrices

  export BLOSUM62_Pi, BLOSUM62_Pij,

  SequenceWeights, Pseudocount, AdditiveSmoothing,
  ResidueContingencyTables, ResidueCount, ResidueProbability,
  nresidues, update!, apply_pseudocount!, count!, normalize!,
  blosum_pseudofrequencies!, apply_pseudofrequencies!, probabilities,
  delete_dimensions!, delete_dimensions,

  AbstractMeasure, SymmetricMeasure, Entropy,
  MutualInformation, MutualInformationOverEntropy,
  estimate, estimate_on_marginal,
  GapUnionPercentage, GapIntersectionPercentage,

  estimateincolumns, estimateinsequences,

  APC!,

  buslje09, BLMI, pairwisegapfraction,

  # Formats from MIToS.MSA:
  Raw, Stockholm, FASTA

  include("BLOSUM62.jl")
  include("Probabilities.jl")
  include("InformationMeasures.jl")
  include("Iterations.jl")
  include("Corrections.jl")
  include("CorrectedMutualInformation.jl")

end
