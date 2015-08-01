using Base.Test
using MIToS.Information
using MIToS.Clustering
using MIToS.MSA

## Test Counts

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

println("### Test count! whit Clusters ###")
false_clusters = Clusters(zeros(20),zeros(20),rand(20))
seq = Residue(MSA._to_char)
@test_throws BoundsError count!(zeros(ResidueCount{Float64, 1, true}), false_clusters, seq)
@test count!(zeros(ResidueCount{Float64, 1, false}), false_clusters, seq).counts == getweight(false_clusters)
@test count!(zeros(ResidueCount{Float64, 2, false}), false_clusters, seq, seq).marginals[:,1] == getweight(false_clusters)
@test count!(zeros(ResidueCount{Float64, 3, false}), false_clusters, seq, seq, seq).marginals[:,1] == getweight(false_clusters)

println("### Test count ###")
@test count(false, false_clusters, seq).counts == getweight(false_clusters)
@test count(false, false_clusters, seq, seq).marginals[:,1] == getweight(false_clusters)
@test count(false, false_clusters, seq, seq, seq).marginals[:,1] == getweight(false_clusters)
@test count(true, seq).total == 21
@test count(true, seq, seq).total == 21
@test count(true, seq, seq, seq).total == 21

## Test Probabilities

@test size(ResidueProbability{1, false}().probabilities) == (20,)
@test size(ResidueProbability{1, true}().probabilities) == (21,)

@test size(ResidueProbability{2, false}().probabilities) == (20,20)
@test size(ResidueProbability{2, true}().probabilities) == (21,21)

@test size(ResidueProbability{3, false}().probabilities) == (20,20,20)
@test size(ResidueProbability{3, true}().probabilities) == (21,21,21)

for ndim = 1:4
	@test size(ResidueProbability{ndim, true}().marginals) == (21, ndim)
	@test size(ResidueProbability{ndim, false}().marginals) == (20, ndim)
end

for ndim = 1:4, usegap = Bool[true, false]
	println("### N: $(ndim) UseGap: $(usegap) ###")
  nres = usegap ? 21 : 20

  println("Test zeros for ResidueProbability{$(ndim), $(usegap)}")
  N = zeros(ResidueProbability{ndim, usegap})
  @test sum(N.probabilities) == zero(Float64)

  println("Test iteration interface for ResidueProbability{$(ndim), $(usegap)}")
  @test collect(N) == zeros(Float64, nres^ndim)

  println("Test size for ResidueProbability{$(ndim), $(usegap)}")
  @test size(N) == size(N.probabilities)

  println("Test indexing with i for ResidueProbability{$(ndim), $(usegap)}")
  @test N[1] == zero(Float64)
  @test N[Int(Residue('R'))] == zero(Float64)

  println("Test setindex! with i for ResidueProbability{$(ndim), $(usegap)}")
  N[3] = 1.0
  @test N[3] == one(Float64)

  println("Test update! for ResidueProbability{$(ndim), $(usegap)}")
  N = zeros(ResidueProbability{ndim, usegap})
  N[1] = 0.5
  update!(N)
  @test N.marginals[1,1] == 0.5

  println("Test normalize! for ResidueProbability{$(ndim), $(usegap)}")
  normalize!(N)
  @test N[1] == 1.0
  @test N.marginals[1,1] == 1.0
  @test sum(N) == 1.0

  println("Test fill! with ResidueCount{Int, $(ndim), $(usegap)} for ResidueProbability{$(ndim), $(usegap)}")
  C = fill!(zeros(ResidueCount{Int, ndim, usegap}), AdditiveSmoothing(1))
  P = fill!(N, C)
  @test P[1] == 1.0/length(C.counts)
  @test P.marginals[1,1] == 1.0 / nres
  @test_approx_eq sum(P) 1.0
end

for usegap = Bool[true, false]
	println("### UseGap: $(usegap) ###")
  nres = usegap ? 21 : 20

	N = zeros(ResidueProbability{2, usegap})

  println("Test indexing with ij for ResidueProbability{2, $(usegap)}")
  @test N[1,1] == zero(Float64)
  @test N[Int(Residue('R')), Int(Residue('R'))] == zero(Float64)

  println("Test setindex! with ij for ResidueProbability{2, $(usegap)}")
  N[3,3] = 1.0
  @test N[3,3] == one(Float64)
end

for ndim = 1:4
	@test size(similar(ResidueProbability{ndim, false}()).marginals) == (20, ndim)
	@test size(similar(ResidueProbability{ndim, true}()).marginals) == (21, ndim)
end

@test size(similar(ResidueProbability{1, false}(),4).marginals) == (20, 4)

# println("### Test probabilities ###")
# @test probabilities(count(false, seq)).marginals[1,1] == 0.05
# @test probabilities(count(false, seq, seq)).marginals[1,1] == 0.05
# @test probabilities(count(false, seq, seq, seq)).marginals[1,1] == 0.05
# @test probabilities(false, seq).marginals[1,1] == 0.05
# @test probabilities(false, seq, seq).marginals[1,1] == 0.05
# @test probabilities(false, seq, seq, seq).marginals[1,1] == 0.05
# @test sum(probabilities(false, seq)) == 1.0
# @test sum(probabilities(false, seq, seq)) == 1.0
# @test sum(probabilities(false, seq, seq, seq)) == 1.0
