module Information

  using MIToS.MSA
  using MIToS.Clustering

  export BLOSUM62_Pi, BLOSUM62_Pij,

  Pseudocount, AdditiveSmoothing,
	ResidueContingencyTables, ResidueCount, ResidueProbability,
	nresidues, update!, apply_pseudocount!, count!, normalize!,
  blosum_pseudofrequencies!, apply_pseudofrequencies!, probabilities,
  delete_dimensions!, delete_dimensions,

  InformationMeasure, Entropy, MutualInformation,
  estimate, estimate_on_marginal

  #Fixed, Pseudofrequencies, ResidueProbabilities, ResiduePairProbabilities

  include("BLOSUM62.jl")
  include("Probabilities.jl")
  include("InformationMeasures.jl")

end
