module Information

  using MIToS.Utils
  using MIToS.MSA
  using MIToS.Clustering
  using PairwiseListMatrices

  export BLOSUM62_Pi, BLOSUM62_Pij,

  SequenceWeights, Pseudocount, AdditiveSmoothing,
  ResidueContingencyTables, ResidueCount, ResidueProbability,
  nresidues, update!, apply_pseudocount!, count!, normalize!,
  blosum_pseudofrequencies!, apply_pseudofrequencies!, probabilities,
  delete_dimensions!, delete_dimensions,

  InformationMeasure, SymmetricMeasure, Entropy,
  MutualInformation, MutualInformationOverEntropy,
  estimate, estimate_on_marginal,
  GapUnionPercentage, GapIntersectionPercentage,

  estimateincolumns, estimateinsequences,

  APC!,

  buslje09, BLMI, pairwisegappercentage,

  # Formats from MIToS.MSA:
  Raw, Stockholm, FASTA

  include("BLOSUM62.jl")
  include("Probabilities.jl")
  include("InformationMeasures.jl")
  include("Iterations.jl")
  include("Corrections.jl")
  include("CorrectedMutualInformation.jl")

end
