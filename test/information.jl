using Base.Test
using MIToS.Information
using MIToS.MSA

@test size(ResidueCount{Int, 1, false}().counts) == (20,)
@test size(ResidueCount{Int, 1, true}().counts) == (21,)

@test size(ResidueCount{Int, 2, false}().counts) == (20,20)
@test size(ResidueCount{Int, 2, true}().counts) == (21,21)

@test size(ResidueCount{Int, 3, false}().counts) == (20,20,20)
@test size(ResidueCount{Int, 3, true}().counts) == (21,21,21)

for ndim = 1:4
	@test size(ResidueCount{Int, ndim, true}().marginals) == (21, ndim)
	@test size(ResidueCount{Int, ndim, false}().marginals) == (20, ndim)
end

for ndim = 1:4, usegap = Bool[true, false]
	println("### N: $(ndim) UseGap: $(usegap) ###")
  nres = usegap ? 21 : 20

  println("Test zeros for ResidueCount{Float64, $(ndim), $(usegap)}")
  N = zeros(ResidueCount{Float64, ndim, usegap})
  @test sum(N.counts) == zero(Float64)

  println("Test iteration interface for ResidueCount{Float64, $(ndim), $(usegap)}")
  @test collect(N) == zeros(Int, nres^ndim)

  println("Test size for ResidueCount{Float64, $(ndim), $(usegap)}")
  @test size(N) == size(N.counts)

  println("Test indexing with i for ResidueCount{Float64, $(ndim), $(usegap)}")
  @test N[1] == zero(Float64)
  @test N[Int(Residue('R'))] == zero(Float64)

  println("Test setindex! with i for ResidueCount{Float64, $(ndim), $(usegap)}")
  N[3] = 1
  @test N[3] == one(Float64)

  println("Test update! for ResidueCount{Float64, $(ndim), $(usegap)}")
  fill!(N.counts, 1)
  update!(N)
  @test N.total == length(N.counts)
  @test N.marginals[1] == nres ^ ( ndim - 1 )

  println("Test apply_pseudocount! with Laplace Smoothing for ResidueCount{Float64, $(ndim), $(usegap)}")
  P = apply_pseudocount!(zeros(ResidueCount{Float64, ndim, usegap}), AdditiveSmoothing(1.0))
  @test P[1,1] == 1.0
  @test P.total == length(P.counts)
  @test P.marginals[1] == nres ^ ( ndim - 1 )

  println("Test fill! with AdditiveSmoothing(0.05) for ResidueCount{Float64, $(ndim), $(usegap)}")
  fill!(P, AdditiveSmoothing(0.05))
  @test P[1,1] == 0.05
  @test P.total == length(P.counts) * 0.05
  @test P.marginals[1] == 0.05 * ( nres ^ ( ndim - 1 ) )
end

for usegap = Bool[true, false]
	println("### UseGap: $(usegap) ###")
  nres = usegap ? 21 : 20
  aas = Residue(MSA._to_char)

	N = zeros(ResidueCount{Float64, 2, usegap})

  println("Test indexing with ij for ResidueCount{Float64, 2, $(usegap)}")
  @test N[1,1] == zero(Float64)
  @test N[Int(Residue('R')), Int(Residue('R'))] == zero(Float64)

  println("Test setindex! with ij for ResidueCount{Float64, 2, $(usegap)}")
  N[3,3] = 1
  @test N[3,3] == one(Float64)

  println("Test count! for ResidueCount{Float64, 2, $(usegap)} and ResidueCount{Float64, 1, $(usegap)}")
	C = zeros(ResidueCount{Float64, 2, usegap})
  count!(C, aas, aas)
  @test C.counts == eye(nres)
  @test (count!(zeros(ResidueCount{Float64, 1, usegap}), aas)).total == nres

  println("Test count! for ResidueCount{Float64, 5, $(usegap)} and  ResidueCount{Float64, 5, $(usegap)}}")
  @test (count!(zeros(ResidueCount{Float64, 5, usegap}), aas, aas, aas, aas, aas)).total == nres

  println("Test count! for ResidueCount{Float64, 4, $(usegap)} and  ResidueCount{Float64, 4, $(usegap)}}")
  C4 = count!(zeros(ResidueCount{Float64, 4, usegap}), aas, aas, aas, aas)
  @test C4[20,20,20,20] == 1.0
  @test C4[20,20,20,19] == 0.0
  # Also test sum
  @test sum(C4) == nres
end

for ndim = 1:4
	@test size(similar(ResidueCount{Int, ndim, false}()).marginals) == (20, ndim)
	@test size(similar(ResidueCount{Int, ndim, true}()).marginals) == (21, ndim)
end

@test size(similar(ResidueCount{Int, 1, false}(),4).marginals) == (20, 4)
@test isa(similar(ResidueCount{Int, 1, false}(),Float64).total, Float64)

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
