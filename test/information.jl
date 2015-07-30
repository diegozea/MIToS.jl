using Base.Test
using MIToS.Information
using MIToS.MSA

@test size(ResidueCount{Int, 1, false}().counts) == (20,)
@test size(ResidueCount{Int, 1, true}().counts) == (21,)

@test size(ResidueCount{Int, 2, false}().counts) == (20,20)
@test size(ResidueCount{Int, 2, true}().counts) == (21,21)

@test size(ResidueCount{Int, 3, false}().counts) == (20,20,20)
@test size(ResidueCount{Int, 3, true}().counts) == (21,21,21)

@test size(ResidueCount{Int, 1, true}().marginals) == (1,21)
@test size(ResidueCount{Int, 1, false}().marginals) == (1,20)

@test size(ResidueCount{Int, 2, true}().marginals) == (2,21)
@test size(ResidueCount{Int, 2, false}().marginals) == (2,20)

@test size(ResidueCount{Int, 3, true}().marginals) == (3,21)
@test size(ResidueCount{Int, 3, false}().marginals) == (3,20)

for ndim = 1:3, usegap = Bool[true, false]
	println("### N: $(ndim) UseGap: $(usegap) ###")
  println("Test zeros for ResidueCount{Float64, $(ndim), $(usegap)}")
  N = zeros(ResidueCount{Float64, ndim, usegap})
  @test sum(N.counts) == zero(Float64)
  println("Test iteration interface for ResidueCount{Float64, $(ndim), $(usegap)}")
  n = usegap ? 21 : 20
  @test collect(N) == zeros(Int, n^ndim)
  println("Test size for ResidueCount{Float64, $(ndim), $(usegap)}")
  @test size(N) == size(N.counts)
  println("Test indexing with i for ResidueCount{Float64, $(ndim), $(usegap)}")
  @test N[1] == zero(Float64)
  @test N[Int(Residue('R'))] == zero(Float64)
end

for usegap = Bool[true, false]
	println("### UseGap: $(usegap) ###")
	N = zeros(ResidueCount{Float64, 2, usegap})
  println("Test indexing with ij for ResidueCount{Float64, 2, $(usegap)}")
  @test N[1,1] == zero(Float64)
  @test N[Int(Residue('R')), Int(Residue('R'))] == zero(Float64)
end



# # Copy & deepcopy
# cnone = copy(none)
# dcnone = deepcopy(none)
# 
# none.Gab[1,1] = 1.0
# 
# @test cnone.α == zero(Float64)
# @test cnone.β == zero(Float64)
# @test cnone.Gab[1,1] == zero(Float64)
# 
# @test dcnone.α == zero(Float64)
# @test dcnone.β == zero(Float64)
# @test dcnone.Gab[1,1] == zero(Float64)
# 
# ## ResidueProbabilities ##
# 
# none = zeros(ResidueProbabilities)
# @test none.Pa == zeros(Float64, 20)
# 
# # Copy & deepcopy
# cnone = copy(none)
# dcnone = deepcopy(none)
# 
# none.Pa[1] = 1.0
# 
# @test cnone.Pa[1] == zero(Float64)
# @test dcnone.Pa[1] == zero(Float64)
# 
# ## ResidueProbabilities ##
# 
# none =  zeros(ResiduePairProbabilities)
# @test none.Pab == zeros(Float64, (20,20))
# @test none.Pa == zeros(Float64, 20)
# @test none.Pb == zeros(Float64, 20)
# 
# # Copy & deepcopy
# cnone = copy(none)
# dcnone = deepcopy(none)
# 
# none.Pab[1,1] = 1.0
# none.Pa[1] = 1.0
# none.Pb[1] = 1.0
# 
# @test cnone.Pab[1,1] == zero(Float64)
# @test cnone.Pa[1] == zero(Float64)
# @test cnone.Pb[1] == zero(Float64)
# 
# @test dcnone.Pab[1,1] == zero(Float64)
# @test dcnone.Pa[1] == zero(Float64)
# @test dcnone.Pb[1] == zero(Float64)
# 
