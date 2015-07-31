module Information

  using MIToS.MSA
  using MIToS.Clustering

  export BLOSUM62_Pa, BLOSUM62_Pab,

  Pseudocount, AdditiveSmoothing, ResidueCount, 
	update!, apply_pseudocount!

  #Fixed, Pseudofrequencies, ResidueProbabilities, ResiduePairProbabilities

  include("BLOSUM62.jl")
  include("Probabilities.jl")

end
