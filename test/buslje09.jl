using Base.Test
using MIToS.Information
using MIToS.Utils
using MIToS.MSA

const Gaoetal2011 = joinpath(pwd(), "data", "Gaoetal2011.fasta")

gao11_buslje09(measure) = joinpath(pwd(), "data", string("data_Gaoetal2011_soft_Busljeetal2009_measure_", measure, ".txt"))

## Column numbers for the output of Buslje et. al. 2009
const SCORE = 9
const ZSCORE = 12

print("""

Simple examples
===============
""")

let aln = Residue[ 'A' 'A'
                   'A' 'R' ]

  estimateincolumns(aln, Float64, ResidueProbability{Float64, 2, false}, MutualInformation{Float64}(), AdditiveSmoothing(0.05))

  print("""
  MI
  """)
  Pij = zeros(Float64, 20, 20);
  N = 2 # There are 2 sequences
  Pij[1,1] = (1/N) # A A
  Pij[1,2] = (1/N) # A R
  @test_approx_eq sum(Pij) 1.0
  Pi = squeeze(Base.sum(Pij,2),2)
  @test Pi[1] == 1.0 # Is always A
  Pj = squeeze(Base.sum(Pij,1),1)
  @test Pj[1] == 0.5 # A
  @test Pj[2] == 0.5 # R
  total = 0.0
  for i in 1:20, j in 1:20
    if Pij[i,j] != 0.0 && Pi[i] != 0.0 && Pj[j] != 0.0
      total += (Pij[i,j] * log(Pij[i,j]/(Pi[i]*Pj[j]))) # 0.5 * log(0.5/(0.5 * 1.0)) == 0.0
    end
  end
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
  cal = [ 0.0  total
        total 0.0  ]
  mi = estimateincolumns(aln, Float64, ResidueProbability{Float64, 2, false}, MutualInformation{Float64}(), AdditiveSmoothing(0.05))
  @test_approx_eq mi[1,2] total

  print("""

  Test APC!
  """)
  mimean = mean(cal) # MI.. == MI.j == MIi. == Mean
  # APC = (mimean * mimean) / mimean == mimean
  APC!(mi)
  @test_approx_eq mi[1,2] mimean

  print("""

  Z-score
  """)
  # 2 Possibilities:  A A   A A
  #                   A R   R A
  # The same MI APC matrix
  mi, zscore = buslje09(aln, lambda=0.05, clustering=false, apc=true)
  println(zscore)

end


print("""

Results from Buslje et. al 2009
===============================
""")

let data = readdlm(gao11_buslje09("MI")); results = buslje09(Gaoetal2011, FASTA, lambda=0.0, clustering=false, apc=false)
  @test_approx_eq_eps maximum(abs(convert(Vector{Float64}, data[:, SCORE]) .- matrix2list(results[2]))) 0.0 0.000001
#  @test cor(convert(Vector{Float64}, data[:, ZSCORE]), matrix2list(results[3])) > 0.90

  println(cor(convert(Vector{Float64}, data[:, ZSCORE]), matrix2list(results[3])))
end



# let data = readdlm(gao11_buslje09("MI_clustering")); results = buslje09(Gaoetal2011, FASTA, lambda=0.0, clustering=true, apc=false)
#   @test_approx_eq_eps maximum(abs(convert(Vector{Float64}, data[:, SCORE]) .- matrix2list(results[2]))) 0.0 0.000001
#   @test cor(convert(Vector{Float64}, data[:, ZSCORE]), matrix2list(results[3])) > 0.90

#   println(cor(convert(Vector{Float64}, data[:, ZSCORE]), matrix2list(results[3])))
# end

# let data = readdlm(gao11_buslje09("MI_APC")); results = buslje09(Gaoetal2011, FASTA, lambda=0.0, clustering=false, apc=true)
#   #@test_approx_eq_eps maximum(abs(convert(Vector{Float64}, data[:, SCORE]) .- matrix2list(results[2]))) 0.0 0.000001
#   @test cor(convert(Vector{Float64}, data[:, ZSCORE]), matrix2list(results[3])) > 0.90

#   println(cor(convert(Vector{Float64}, data[:, ZSCORE]), matrix2list(results[3])))
# end

