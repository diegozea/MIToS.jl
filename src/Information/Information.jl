module Information

  using MIToS.MSA
  using MIToS.Clustering

  export BLOSUM62_Pa, BLOSUM62_Pab,

  Pseudocount, AdditiveSmoothing, ResidueCount, ResidueProbability,
	update!, apply_pseudocount!, count!


  #Fixed, Pseudofrequencies, ResidueProbabilities, ResiduePairProbabilities

  include("BLOSUM62.jl")
  include("Probabilities.jl")

end
