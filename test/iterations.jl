using Base.Test

using MIToS.MSA
using MIToS.Information

let aln = read("./data/Gaoetal2011.fasta", FASTA), result = Float64[ 0     0     0     0     0     0
                                                                     0     0     0     0     0     0
                                                                     0     0     0     1     1     0.296
                                                                     0     0     1     0     1     0.296
                                                                     0     0     1     1     0     0.296
                                                                     0     0     0.296 0.296 0.296 0 ]

@test_approx_eq_eps estimateincolumns(aln, ResidueCount{Int, 2, false}, MutualInformationOverEntropy{Float64}()) result 0.0001
@test_approx_eq_eps estimateincolumns(aln, Int, ResidueProbability{Float64, 2, false}, MutualInformationOverEntropy{Float64}()) result 0.0001
