# using Base.Test
# using MIToS.MSA
# using MIToS.Information

print("""

Estimation of measures on a MSA
===============================
""")

print("""
This is the example of MI(X, Y)/H(X, Y) from:
Gao, H., Dou, Y., Yang, J., & Wang, J. (2011). New methods to measure residues coevolution in proteins. BMC bioinformatics, 12(1), 206.
""")
let aln = read(joinpath(pwd(), "data", "Gaoetal2011.fasta"), FASTA), result = Float64[ 0     0     0     0     0     0
                                                                     0     0     0     0     0     0
                                                                     0     0     0     1     1     0.296
                                                                     0     0     1     0     1     0.296
                                                                     0     0     1     1     0     0.296
                                                                     0     0     0.296 0.296 0.296 0 ]

  @test_approx_eq_eps triu(estimateincolumns(aln, ResidueCount{Int, 2, false}, MutualInformationOverEntropy{Float64}()), 1) triu(result, 1) 0.0001
  @test_approx_eq_eps triu(estimateincolumns(aln, Int, ResidueProbability{Float64, 2, false}, MutualInformationOverEntropy{Float64}()), 1) triu(result, 1) 0.0001
  @test_approx_eq_eps triu(estimateincolumns(aln, Float64, ResidueProbability{Float64, 2, false}, MutualInformationOverEntropy{Float64}()), 1) triu(result, 1) 0.0001
  @test_approx_eq_eps triu(estimateincolumns(aln, BigFloat, ResidueProbability{BigFloat, 2, false}, MutualInformationOverEntropy{BigFloat}()), 1) triu(result, 1) 0.0001

  @test_approx_eq_eps estimateincolumns(aln, ResidueCount{Int, 2, false}, MutualInformationOverEntropy{Float64}(),
                                        zero(AdditiveSmoothing{Int}), NoClustering(), false, zero(Float64)) result 0.0001
  @test_approx_eq_eps estimateincolumns(aln, Int, ResidueProbability{Float64, 2, false}, MutualInformationOverEntropy{Float64}(),
                                        zero(AdditiveSmoothing{Int}), NoClustering(), false, zero(Float64)) result 0.0001
  @test_approx_eq_eps estimateincolumns(aln, Float64, ResidueProbability{Float64, 2, false}, MutualInformationOverEntropy{Float64}(),
                                        zero(AdditiveSmoothing{Float64}), NoClustering(), false, zero(Float64)) result 0.0001
  @test_approx_eq_eps estimateincolumns(aln, BigFloat, ResidueProbability{BigFloat, 2, false}, MutualInformationOverEntropy{BigFloat}(),
                                        zero(AdditiveSmoothing{BigFloat}), NoClustering(), false, zero(BigFloat)) result 0.0001
end

