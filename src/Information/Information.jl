module Information

  using MIToS.MSA
  using MIToS.Clustering

  export BLOSUM62_Pi, BLOSUM62_Pij,

  Pseudocount, AdditiveSmoothing, 
	ResidueContingencyTables, ResidueCount, ResidueProbability,
	nresidues, update!, apply_pseudocount!, count!, normalize!,
  blosum_pseudofrequencies!, apply_pseudofrequencies!


  #Fixed, Pseudofrequencies, ResidueProbabilities, ResiduePairProbabilities

  include("BLOSUM62.jl")
  include("Probabilities.jl")

end
