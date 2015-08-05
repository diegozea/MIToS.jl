using Base.Test
using MIToS.Information
using MIToS.Clustering
using MIToS.MSA

const false_clusters = Clusters(zeros(20),zeros(20),rand(20))
const seq = Residue(MSA._to_char)

## Test Counts

@test size(ResidueCount{Int, 1, false}().table) == (20,)
@test size(ResidueCount{Int, 1, true}().table) == (21,)

@test size(ResidueCount{Int, 2, false}().table) == (20,20)
@test size(ResidueCount{Int, 2, true}().table) == (21,21)

@test size(ResidueCount{Int, 3, false}().table) == (20,20,20)
@test size(ResidueCount{Int, 3, true}().table) == (21,21,21)

for ndim = 1:4
	@test size(ResidueCount{Int, ndim, true}().marginals) == (21, ndim)
	@test size(ResidueCount{Int, ndim, false}().marginals) == (20, ndim)
end

for ndim = 1:4, usegap = Bool[true, false]
	println("### N: $(ndim) UseGap: $(usegap) ###")
  nres = usegap ? 21 : 20

	println("Test nresidues for ResidueCount{Float64, $(ndim), $(usegap)}")
	N = zeros(ResidueCount{Float64, ndim, usegap})
	@test nresidues(N) == nres

  println("Test zeros for ResidueCount{Float64, $(ndim), $(usegap)}")
  @test sum(N.table) == zero(Float64)

  println("Test iteration interface for ResidueCount{Float64, $(ndim), $(usegap)}")
  @test collect(N) == zeros(Int, nres^ndim)

  println("Test size for ResidueCount{Float64, $(ndim), $(usegap)}")
  @test size(N) == size(N.table)

  println("Test indexing with i for ResidueCount{Float64, $(ndim), $(usegap)}")
  @test N[1] == zero(Float64)
  @test N[Int(Residue('R'))] == zero(Float64)

  println("Test setindex! with i for ResidueCount{Float64, $(ndim), $(usegap)}")
  N[3] = 1
  @test N[3] == one(Float64)

  println("Test update! for ResidueCount{Float64, $(ndim), $(usegap)}")
  fill!(N.table, 1)
  update!(N)
  @test N.total == length(N.table)
  @test N.marginals[1] == nres ^ ( ndim - 1 )

  println("Test apply_pseudocount! with Laplace Smoothing for ResidueCount{Float64, $(ndim), $(usegap)}")
  P = apply_pseudocount!(zeros(ResidueCount{Float64, ndim, usegap}), AdditiveSmoothing(1.0))
  @test P[1,1] == 1.0
  @test P.total == length(P.table)
  @test P.marginals[1] == nres ^ ( ndim - 1 )

  println("Test fill! with AdditiveSmoothing(0.05) for ResidueCount{Float64, $(ndim), $(usegap)}")
  fill!(P, AdditiveSmoothing(0.05))
  @test P[1,1] == 0.05
  @test P.total == length(P.table) * 0.05
  @test P.marginals[1] == 0.05 * ( nres ^ ( ndim - 1 ) )
end

for usegap = Bool[true, false]
	println("### UseGap: $(usegap) ###")
  nres = usegap ? 21 : 20

	N = zeros(ResidueCount{Float64, 2, usegap})

  println("Test indexing with ij for ResidueCount{Float64, 2, $(usegap)}")
  @test N[1,1] == zero(Float64)
  @test N[Int(Residue('R')), Int(Residue('R'))] == zero(Float64)

  println("Test setindex! with ij for ResidueCount{Float64, 2, $(usegap)}")
  N[3,3] = 1
  @test N[3,3] == one(Float64)

  println("Test count! for ResidueCount{Float64, 2, $(usegap)} and ResidueCount{Float64, 1, $(usegap)}")
	C = zeros(ResidueCount{Float64, 2, usegap})
  count!(C, seq, seq)
  @test C.table == eye(nres)
  @test (count!(zeros(ResidueCount{Float64, 1, usegap}), seq)).total == nres

  println("Test count! for ResidueCount{Float64, 5, $(usegap)} and  ResidueCount{Float64, 5, $(usegap)}}")
  @test (count!(zeros(ResidueCount{Float64, 5, usegap}), seq, seq, seq, seq, seq)).total == nres

  println("Test count! for ResidueCount{Float64, 4, $(usegap)} and  ResidueCount{Float64, 4, $(usegap)}}")
  C4 = count!(zeros(ResidueCount{Float64, 4, usegap}), seq, seq, seq, seq)
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
@test_throws BoundsError count!(zeros(ResidueCount{Float64, 1, true}), false_clusters, seq)
@test count!(zeros(ResidueCount{Float64, 1, false}), false_clusters, seq).table == getweight(false_clusters)
@test count!(zeros(ResidueCount{Float64, 2, false}), false_clusters, seq, seq).marginals[:,1] == getweight(false_clusters)
@test count!(zeros(ResidueCount{Float64, 3, false}), false_clusters, seq, seq, seq).marginals[:,1] == getweight(false_clusters)

println("### Test count ###")
@test count(seq, weight=false_clusters).table == getweight(false_clusters)
@test count(seq, seq, weight=false_clusters).marginals[:,1] == getweight(false_clusters)
@test count(seq, seq, seq, weight=false_clusters).marginals[:,1] == getweight(false_clusters)

@test count(seq, usegap=true).total == 21
@test count(seq, seq, usegap=true).total == 21
@test count(seq, seq, seq, usegap=true).total == 21

@test count(AdditiveSmoothing(1.0), seq, usegap=true).total == 21 + 21
@test count(AdditiveSmoothing(1.0), seq, seq, usegap=true).total == 21 + (21*21)


## Test Probabilities

@test size(ResidueProbability{1, false}().table) == (20,)
@test size(ResidueProbability{1, true}().table) == (21,)

@test size(ResidueProbability{2, false}().table) == (20,20)
@test size(ResidueProbability{2, true}().table) == (21,21)

@test size(ResidueProbability{3, false}().table) == (20,20,20)
@test size(ResidueProbability{3, true}().table) == (21,21,21)

for ndim = 1:4
	@test size(ResidueProbability{ndim, true}().marginals) == (21, ndim)
	@test size(ResidueProbability{ndim, false}().marginals) == (20, ndim)
end

for ndim = 1:4, usegap = Bool[true, false]
	println("### N: $(ndim) UseGap: $(usegap) ###")
  nres = usegap ? 21 : 20

  println("Test zeros for ResidueProbability{$(ndim), $(usegap)}")
  N = zeros(ResidueProbability{ndim, usegap})
  @test sum(N.table) == zero(Float64)

  println("Test iteration interface for ResidueProbability{$(ndim), $(usegap)}")
  @test collect(N) == zeros(Float64, nres^ndim)

  println("Test size for ResidueProbability{$(ndim), $(usegap)}")
  @test size(N) == size(N.table)

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
  @test P[1] == 1.0/length(C.table)
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

println("### Test pseudofrequencies ###")

Pab = fill!(zeros(ResidueProbability{2, false}), count(seq, seq))

Gab = blosum_pseudofrequencies!(zeros(ResidueProbability{2, false}), Pab)
@test_approx_eq sum(Gab) 1.0

@test_approx_eq sum( apply_pseudofrequencies!(copy(Pab), Gab, 1, 1) ) 1.0

Pab_with_pseudofrequency = apply_pseudofrequencies!(copy(Pab), Gab, 10, 0)
@test Pab_with_pseudofrequency == Pab


println("### Test probabilities ###")
@test probabilities(seq, weight=false_clusters).table == getweight(false_clusters)/sum(getweight(false_clusters))
@test probabilities(seq, seq, weight=false_clusters).marginals[:,1] == getweight(false_clusters)/sum(getweight(false_clusters))
@test probabilities(seq, seq, seq, weight=false_clusters).marginals[:,1] == getweight(false_clusters)/sum(getweight(false_clusters))

@test_approx_eq probabilities(seq)[1] 0.05
@test_approx_eq probabilities(seq, seq)[2,2] 0.05
@test_approx_eq probabilities(seq, seq, seq)[3,3,3] 0.05

@test_approx_eq probabilities(AdditiveSmoothing(1.0), seq, usegap=true)[1] 1/21.0
Pab = probabilities(AdditiveSmoothing(1.0), seq, seq, usegap=true)
@test_approx_eq Pab[1,1] 2.0/(21 + 21*21)
@test_approx_eq Pab[1,2] 1.0/(21 + 21*21)
@test_approx_eq sum(Pab) 1.0

@test_approx_eq probabilities(1, 0, seq, seq)[2,2] 0.05

println("### Test delete_dimensions ###")
Pxyz = probabilities(seq, seq, seq);
@test delete_dimensions(Pxyz, 3) == probabilities(seq, seq)
@test delete_dimensions(Pxyz, 3, 2) == probabilities(seq)
Pxy = delete_dimensions(Pxyz, 3);
@test delete_dimensions!(Pxy, Pxyz, 1) == probabilities(seq, seq)
@test_approx_eq sum(Pxy) 1.0
Nxyz = count(seq, seq, seq);
Nxy = delete_dimensions(Nxyz, 3);
@test Nxy == count(seq, seq)
@test delete_dimensions(Nxyz, 3, 2) == count(seq)
@test delete_dimensions!(Nxy, Nxyz, 1) == count(seq, seq)
@test sum(Nxy) == Nxyz.total



