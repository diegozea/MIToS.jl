# using Base.Test
# using MIToS.Information
# using MIToS.Utils
# using MIToS.MSA
# using PairwiseListMatrices

const Gaoetal2011 = joinpath(pwd(), "data", "Gaoetal2011.fasta")

gao11_buslje09(measure) = joinpath(pwd(), "data", string("data_Gaoetal2011_soft_Busljeetal2009_measure_", measure, ".txt"))

## Column numbers for the output of Buslje et. al. 2009
const SCORE = 9
const ZSCORE = 12

const MIToS_SCORE = 2
const MIToS_ZSCORE = 1

print("""

Test APC!
=========
""")

let MI = [ 0 2 4
           2 0 6
           4 6 0 ]

  mean_col = Information._mean_column(MI)
  mean_tot = Information._mean_total(MI)
  @test mean_col == [3., 4., 5.]
  @test mean_tot == 4.

  MIp = APC!(convert(Matrix{Float64}, MI))
  @test_approx_eq MIp [  NaN -1.0  0.25
                        -1.0  NaN  1.00
                        0.25 1.00   NaN ]
end

let MI = PairwiseListMatrix([2, 4, 6])

  mean_col = mean_nodiag(MI, 1)
  mean_tot = mean_nodiag(MI)
  @test mean_col == [3. 4. 5.]
  @test mean_tot == 4.

  MIp = APC!(convert(PairwiseListMatrix{Float64, false}, MI))
  @test_approx_eq MIp [  NaN -1.0  0.25
                        -1.0  NaN  1.00
                        0.25 1.00   NaN ]
end

let MI = PairwiseListMatrix([0, 2, 4, 0, 6, 0], true)

  mean_col = mean_nodiag(MI, 1)
  mean_tot = mean_nodiag(MI)
  @test mean_col == [3. 4. 5.]
  @test mean_tot == 4.

  MIp = APC!(convert(PairwiseListMatrix{Float64, true}, MI))
  @test_approx_eq MIp [  NaN -1.0  0.25
                        -1.0  NaN  1.00
                        0.25 1.00   NaN ]
end

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

    # MI.. == MI.j == MIi. == total
    # APC = (total * total) / total == total
    # MI - APC == total - total == 0.0
    APC!(zerodiagonal)
    @test_approx_eq zerodiagonal[1,2] 0.0
  end

  APC!(mi)
  @test isnan(mi[1,1])
  @test_approx_eq_eps mi[1,2] 0.0 1e-16

  print("""

  Z-score
  """)
  # 2 Possibilities:  A A   A A
  #                   A R   R A
  # Is almost the same MI, Z score should be 0.0
  results = buslje09(aln, lambda=0.05, clustering=false, apc=false)
  @test_approx_eq results[MIToS_ZSCORE][1,2] 0.0

  @test aln == Residue[ 'A' 'A'
                        'A' 'R' ]

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

    # MI.. == MI.j == MIi. == total
    # APC = (total * total) / total == total
    # MI - APC == total - total == 0.0
    APC!(zerodiagonal)
    @test_approx_eq zerodiagonal[1,2] 0.0
  end

  APC!(mi)
  @test isnan(mi[1,1])
  @test_approx_eq mi[1,2] 0.0

  print("""

  Z-score
  """)
  # 4 Possibilities:  R A   R A   A R   A R
  #                   A R   R A   A R   R A
  mi = estimateincolumns(aln, Float64, ResidueProbability{Float64, 2, false}, MutualInformation{Float64}(), AdditiveSmoothing(0.05))
  other_posib = Residue[ 'A' 'R'
                         'A' 'R' ]
  other_mi = estimateincolumns(other_posib, Float64, ResidueProbability{Float64, 2, false}, MutualInformation{Float64}(), AdditiveSmoothing(0.05))
  # The mean should be:
  r_mean = 0.5 * ( mi[1,2] + other_mi[1,2] )
  # And the std should be similar to:
  r_std = 0.5 * ( sqrt((mi[1,2]-r_mean)^2) + sqrt((other_mi[1,2]-r_mean)^2) )

  results = buslje09(aln, lambda=0.05, clustering=false, apc=false, samples=100)
  @test_approx_eq_eps results[MIToS_ZSCORE][1,2] ((mi[1,2] - r_mean) / r_std) 0.5

  results = buslje09(aln, lambda=0.05, clustering=false, apc=false, samples=1000)
  @test_approx_eq_eps results[MIToS_ZSCORE][1,2] ((mi[1,2] - r_mean) / r_std) 0.1

  results = buslje09(aln, lambda=0.05, clustering=false, apc=false, samples=10000)
  @test_approx_eq_eps results[MIToS_ZSCORE][1,2] ((mi[1,2] - r_mean) / r_std) 0.05

  @test aln == Residue[ 'R' 'A'
                        'A' 'R' ]

end

print("""

Results from Buslje et. al 2009
===============================
""")

println("""

MI
""")

let data = readdlm(joinpath(pwd(), "data", "data_simple_soft_Busljeetal2009_measure_MI.txt")); results = buslje09(joinpath(pwd(),
           "data", "simple.fasta"), FASTA, lambda=0.0, clustering=false, apc=false)

  @test_approx_eq_eps Float64(data[1, SCORE]) results[MIToS_SCORE][1,2] 1e-6
  @test_approx_eq_eps Float64(data[1, ZSCORE]) results[MIToS_ZSCORE][1,2]  1.5
end

let data = readdlm(gao11_buslje09("MI")); results = buslje09(Gaoetal2011, FASTA, lambda=0.0, clustering=false, apc=false)

  @test_approx_eq_eps convert(Vector{Float64}, data[:, SCORE]) matrix2list(results[MIToS_SCORE]) 1e-6
  @test_approx_eq_eps convert(Vector{Float64}, data[:, ZSCORE]) matrix2list(results[MIToS_ZSCORE]) 1.5

  println("Pearson for Z-score: ", cor(convert(Vector{Float64}, data[:, ZSCORE]), matrix2list(results[MIToS_ZSCORE])))
end

@test_approx_eq_eps buslje09(Gaoetal2011, FASTA, lambda=0.05, clustering=false, apc=false)[MIToS_SCORE][1,2] 0.33051006116310444 1e-14

println("""

MI + clustering
""")

let data = readdlm(gao11_buslje09("MI_clustering")); results = buslje09(Gaoetal2011, FASTA, lambda=0.0, clustering=true, apc=false)

  @test_approx_eq_eps convert(Vector{Float64}, data[:, SCORE]) matrix2list(results[MIToS_SCORE]) 1e-6
  @test_approx_eq_eps convert(Vector{Float64}, data[:, ZSCORE]) matrix2list(results[MIToS_ZSCORE]) 1.5

  println("Pearson for Z-score: ", cor(convert(Vector{Float64}, data[:, ZSCORE]), matrix2list(results[MIToS_ZSCORE])))
end

println("""

MIp
""")

let data = readdlm(gao11_buslje09("MI_APC")); results = buslje09(Gaoetal2011, FASTA, lambda=0.0, clustering=false, apc=true)

  @test_approx_eq_eps convert(Vector{Float64}, data[:, SCORE]) matrix2list(results[MIToS_SCORE]) 1e-6
  @test_approx_eq_eps convert(Vector{Float64}, data[:, ZSCORE]) matrix2list(results[MIToS_ZSCORE]) 1.5

  @test_approx_eq_eps results[MIToS_SCORE][5,6] 0.018484 0.000001

  println("Pearson for MIp: ", cor(convert(Vector{Float64}, data[:, SCORE]), matrix2list(results[MIToS_SCORE])))
  println("Pearson for Z-score: ", cor(convert(Vector{Float64}, data[:, ZSCORE]), matrix2list(results[MIToS_ZSCORE])))
end

println("""

MIp + clustering
""")

let data = readdlm(gao11_buslje09("MI_APC_clustering")); results = buslje09(Gaoetal2011, FASTA, lambda=0.0, clustering=true, apc=true)

  @test_approx_eq_eps convert(Vector{Float64}, data[:, SCORE]) matrix2list(results[MIToS_SCORE]) 1e-6
  @test_approx_eq_eps convert(Vector{Float64}, data[:, ZSCORE]) matrix2list(results[MIToS_ZSCORE]) 1.5

  @test_approx_eq_eps results[MIToS_SCORE][5,6] 0.018484 0.000001

  println("Pearson for MIp: ", cor(convert(Vector{Float64}, data[:, SCORE]), matrix2list(results[MIToS_SCORE])))
  println("Pearson for Z-score: ", cor(convert(Vector{Float64}, data[:, ZSCORE]), matrix2list(results[MIToS_ZSCORE])))
end

# TO DO: Test labels!

print("""

Test for BLMI
=============
""")

let file = joinpath(pwd(), "data", "simple.fasta"),
    busl = buslje09(file, FASTA),
    blmi = BLMI(file, FASTA)

  @test_approx_eq busl[1] blmi[1]
  @test_approx_eq busl[2] blmi[2]

end


let msa = read(Gaoetal2011, FASTA),
    busl = buslje09(Gaoetal2011, FASTA, lambda=0.0, samples=0),
    blmi = BLMI(msa, lambda=0.0, beta=0.0, samples=5) # BLMI should be equal to Buslje09 if beta is zero

  @test_approx_eq busl[2] blmi[2] # MIapc

  @test msa == read(Gaoetal2011, FASTA)

end

let busl = buslje09(Gaoetal2011, FASTA, lambda=0.5, samples=0),
    blmi = BLMI(Gaoetal2011, FASTA, lambda=0.5, beta=0.0, samples=0) # BLMI should be equal to Buslje09 if beta is zero

  @test_approx_eq busl[2] blmi[2] # MIapc
end

print("""

Test for Pairwise Gap Percentage
================================
""")

let file = joinpath(pwd(), "data", "simple.fasta"),
    mat = [ 0. 0.
            0. 0. ]

  (gu, gi) = pairwisegappercentage(file, FASTA)

  @test gu == mat
  @test gi == mat
end

let file = joinpath(pwd(), "data", "gaps.txt")

  gu, gi = pairwisegappercentage(file, Raw)
  cl = hobohmI(read(file, Raw), 0.62)
  ncl = getnclusters(cl)

  @test_approx_eq gu[1, 1] 0.0
  @test_approx_eq gi[1, 1] 0.0

  @test_approx_eq gu[1, 2] getweight(cl, 10)/ncl
  @test_approx_eq gi[1, 2] 0.0

  @test_approx_eq gu[10, 9] (ncl - getweight(cl, 1))/ncl
  @test_approx_eq gi[10, 9] (ncl - getweight(cl, 1) - getweight(cl, 2))/ncl

  @test_approx_eq gu[10, 10] (ncl - getweight(cl, 1))/ncl
  @test_approx_eq gu[10, 10] (ncl - getweight(cl, 1))/ncl
end
