using Base.Test
using MIToS.Information
using MIToS.Utils
using MIToS.MSA

const Gaoetal2011 = joinpath(pwd(), "data", "Gaoetal2011.fasta")

gao11_buslje09(measure) = joinpath(pwd(), "data", string("data_Gaoetal2011_soft_Busljeetal2009_measure_", measure, ".txt"))

## Column numbers for the output of Buslje et. al. 2009
const SCORE = 9
const ZSCORE = 12

const MIToS_SCORE = 2
const MIToS_ZSCORE = 5

print("""

Simple examples
===============
""")

let aln = Residue[ 'A' 'A'
                   'A' 'R' ]

  print("""
  MI
  """)
  # Fill the Pij matrix of the example
  Pij = zeros(Float64, 20, 20);
  N = 2 # There are 2 sequences
  Pij[1,1] = (1/N) # A A
  Pij[1,2] = (1/N) # A R
  @test_approx_eq sum(Pij) 1.0 # Is the matrix correct?
  # Fill the marginals
  Pi = squeeze(Base.sum(Pij,2),2)
  @test Pi[1] == 1.0 # Is always A
  Pj = squeeze(Base.sum(Pij,1),1)
  @test Pj[1] == 0.5 # A
  @test Pj[2] == 0.5 # R
  # Start to sum with 0.0
  total = 0.0
  for i in 1:20, j in 1:20
    if Pij[i,j] != 0.0 && Pi[i] != 0.0 && Pj[j] != 0.0
      total += (Pij[i,j] * log(Pij[i,j]/(Pi[i]*Pj[j]))) # 0.5 * log(0.5/(0.5 * 1.0)) == 0.0
    end
  end
  # Compare with MIToS result
  mi = estimateincolumns(aln, Float64, ResidueProbability{Float64, 2, false}, MutualInformation{Float64}())
  @test mi[1,2] == total

  print("""

  MI: Using pseudocount (0.05)
  ----------------------------
  """)
  Pij = zeros(Float64, 20, 20);
  N = (400 * 0.05) + 2 # 0.05 is the pseudocount and there are 2 sequences
  fill!(Pij, 0.05/N);
  Pij[1,1] =(1.05/N) # A A
  Pij[1,2] =(1.05/N) # A R
  @test_approx_eq sum(Pij) 1.0
  Pi = squeeze(Base.sum(Pij,2),2)
  Pj = squeeze(Base.sum(Pij,1),1)
  total = 0.0
  for i in 1:20, j in 1:20
    total += (Pij[i,j] * log(Pij[i,j]/(Pi[i]*Pj[j])))
  end
  mi = estimateincolumns(aln, Float64, ResidueProbability{Float64, 2, false}, MutualInformation{Float64}(), AdditiveSmoothing(0.05))
  @test_approx_eq mi[1,2] total

  print("""

  Test APC!
  """)
  let zerodiagonal = Float64[ 0.0 total
                             total 0.0 ]

    mimean = mean(zerodiagonal) # MI.. == MI.j == MIi. == Mean
    # APC = (mimean * mimean) / mimean == mimean
    APC!(zerodiagonal)
    @test_approx_eq zerodiagonal[1,2] mimean
  end

  print("""

  Z-score
  """)
  # 2 Possibilities:  A A   A A
  #                   A R   R A
  # Is almost the same MI, Z score should be 0.0
  results = buslje09(aln, lambda=0.05, clustering=false, apc=true)
  @test results[MIToS_ZSCORE] == [ 0.0 0.0
                                   0.0 0.0 ]

end

let aln = Residue[ 'R' 'A'
                   'A' 'R' ]

  print("""
  MI
  """)
  Pij = zeros(Float64, 20, 20);
  N = 2 # There are 2 sequences
  Pij[2,1] = (1/N) # R A
  Pij[1,2] = (1/N) # A R
  @test_approx_eq sum(Pij) 1.0
  Pi = squeeze(Base.sum(Pij,2),2)
  @test Pi[2] == 0.5 # R
  @test Pi[1] == 0.5 # A
  Pj = squeeze(Base.sum(Pij,1),1)
  @test Pj[1] == 0.5 # A
  @test Pj[2] == 0.5 # R
  total = 0.0
  for i in 1:20, j in 1:20
    if Pij[i,j] != 0.0 && Pi[i] != 0.0 && Pj[j] != 0.0
      total += (Pij[i,j] * log(Pij[i,j]/(Pi[i]*Pj[j])))
    end
  end
  @test total == 0.5 * log(2) + 0.5 * log(2) # 0.5/(0.5*0.5) == 2
  mi = estimateincolumns(aln, Float64, ResidueProbability{Float64, 2, false}, MutualInformation{Float64}())
  @test mi[1,2] == total

  print("""

  MI: Using pseudocount (0.05)
  ----------------------------
  """)
  Pij = zeros(Float64, 20, 20);
  N = (400 * 0.05) + 2 # 0.05 is the pseudocount and there are 2 sequences
  fill!(Pij, 0.05/N);
  Pij[2,1] =(1.05/N) # R A
  Pij[1,2] =(1.05/N) # A R
  @test_approx_eq sum(Pij) 1.0
  Pi = squeeze(Base.sum(Pij,2),2)
  Pj = squeeze(Base.sum(Pij,1),1)
  total = 0.0
  for i in 1:20, j in 1:20
    total += (Pij[i,j] * log(Pij[i,j]/(Pi[i]*Pj[j])))
  end
  mi = estimateincolumns(aln, Float64, ResidueProbability{Float64, 2, false}, MutualInformation{Float64}(), AdditiveSmoothing(0.05))
  @test_approx_eq mi[1,2] total

  print("""

  Test APC!
  """)
  let zerodiagonal = Float64[ 0.0 total
                             total 0.0 ]

    mimean = mean(zerodiagonal) # MI.. == MI.j == MIi. == Mean
    # APC = (mimean * mimean) / mimean == mimean
    APC!(zerodiagonal)
    @test_approx_eq zerodiagonal[1,2] mimean
  end

  APC!(mi)

  print("""

  Z-score
  """)
  # 4 Possibilities:  R A   R A   A R   A R
  #                   A R   R A   A R   R A
  other_posib = Residue[ 'A' 'R'
                         'A' 'R' ]
  other_mi = estimateincolumns(other_posib, Float64, ResidueProbability{Float64, 2, false}, MutualInformation{Float64}(), AdditiveSmoothing(0.05))
  APC!(other_mi)
  # The mean should be:
  r_mean = 0.5 * ( mi[1,2] + other_mi[1,2] )
  @test_approx_eq_eps mean(Float64[[mi[1,2], other_mi[1,2]][rand(1:2)] for i in 1:100]) r_mean 1e-15
  # And the std should be similar to:
  r_std = 0.5 * ( sqrt((mi[1,2]-r_mean)^2) + sqrt((other_mi[1,2]-r_mean)^2) )
  @test_approx_eq_eps std(Float64[[mi[1,2], other_mi[1,2]][rand(1:2)] for i in 1:100]) r_std 1e-15
  @test isapprox(r_std + 1.0, 1.0) # This is approx 0.0, so...
  results = buslje09(aln, lambda=0.05, clustering=false, apc=true)
  @test isnan(results[MIToS_ZSCORE][1,2]) # ...the Z score is NaN

end

print("""

Results from Buslje et. al 2009
===============================
""")

let data = readdlm(joinpath(pwd(), "data", "data_simple_soft_Busljeetal2009_measure_MI.txt")); results = buslje09(joinpath(pwd(),
           "data", "simple.fasta"), FASTA, lambda=0.0, clustering=false, apc=false)

  @test_approx_eq_eps Float64(data[1, SCORE]) results[MIToS_SCORE][1,2]  0.000001
  @test_approx_eq_eps Float64(data[1, ZSCORE]) results[MIToS_ZSCORE][1,2]  0.5

  println("Pearson for Z-score: ", cor(convert(Vector{Float64}, data[:, ZSCORE]), matrix2list(results[MIToS_ZSCORE])))
end

let data = readdlm(gao11_buslje09("MI")); results = buslje09(Gaoetal2011, FASTA, lambda=0.0, clustering=false, apc=false)

  @test_approx_eq_eps convert(Vector{Float64}, data[:, SCORE]) matrix2list(results[MIToS_SCORE]) 1e-6
  @test_approx_eq_eps convert(Vector{Float64}, data[:, ZSCORE]) matrix2list(results[MIToS_ZSCORE]) 1.0

  println("Pearson for Z-score: ", cor(convert(Vector{Float64}, data[:, ZSCORE]), matrix2list(results[MIToS_ZSCORE])))
end

let data = readdlm(gao11_buslje09("MI_clustering")); results = buslje09(Gaoetal2011, FASTA, lambda=0.0, clustering=true, apc=false)

  @test_approx_eq_eps convert(Vector{Float64}, data[:, SCORE]) matrix2list(results[MIToS_SCORE]) 1e-6
  @test_approx_eq_eps convert(Vector{Float64}, data[:, ZSCORE]) matrix2list(results[MIToS_ZSCORE]) 1.0

  println("Pearson for Z-score: ", cor(convert(Vector{Float64}, data[:, ZSCORE]), matrix2list(results[MIToS_ZSCORE])))
 end

# # APC is very different now we are using the diagonal

# let data = readdlm(gao11_buslje09("MI_APC")); results = buslje09(Gaoetal2011, FASTA, lambda=0.0, clustering=false, apc=true)
#   #@test_approx_eq_eps maximum(abs(convert(Vector{Float64}, data[:, SCORE]) .- matrix2list(results[2]))) 0.0 0.000001
#   @test cor(convert(Vector{Float64}, data[:, ZSCORE]), matrix2list(results[3])) > 0.90

#   println(cor(convert(Vector{Float64}, data[:, ZSCORE]), matrix2list(results[3])))
# end

