# using Base.Test
# using MIToS.Information
# using MIToS.Clustering
# using MIToS.MSA

let false_clusters = Clusters(zeros(20),zeros(20),rand(20)),
  seq = res"ARNDCQEGHILKMFPSTWYV-"

  print("""

Tests for ResidueContingencyTables
==================================
""")

  print("""

Tests for ResidueCount
----------------------
""")

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

  print("""

Iterate for N & UseGap:
""")
  for ndim = 1:4, usegap = Bool[true, false]
    println("# Iteration = N: $(ndim) UseGap: $(usegap)")
    nres = usegap ? 21 : 20

    println("Test nresidues for ResidueCount{Float64, $(ndim), $(usegap)}")
    N = zeros(ResidueCount{Float64, ndim, usegap})
    @test nresidues(N) == nres

    println("Test zeros for ResidueCount{Float64, $(ndim), $(usegap)}")
    @test sum(N.table) == zero(BigFloat)

    println("Test iteration interface for ResidueCount{Float64, $(ndim), $(usegap)}")
    @test collect(N) == zeros(Int, nres^ndim)

    println("Test size for ResidueCount{Float64, $(ndim), $(usegap)}")
    @test size(N) == size(N.table)

    println("Test indexing with i for ResidueCount{Float64, $(ndim), $(usegap)}")
    @test N[1] == zero(BigFloat)
    @test N[Int(Residue('R'))] == zero(BigFloat)

    println("Test setindex! with i for ResidueCount{Float64, $(ndim), $(usegap)}")
    N[3] = 1
    @test N[3] == one(BigFloat)

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
  print("""
Iteration end

""")

  print("""

Iterate for UseGap in [true, false]:
""")
  for usegap = Bool[true, false]
    println("# Iteration = UseGap: $(usegap)")
    nres = usegap ? 21 : 20

    N = zeros(ResidueCount{Float64, 2, usegap})

    println("Test indexing with ij for ResidueCount{Float64, 2, $(usegap)}")
    @test N[1,1] == zero(BigFloat)
    @test N[Int(Residue('R')), Int(Residue('R'))] == zero(BigFloat)

    println("Test setindex! with ij for ResidueCount{Float64, 2, $(usegap)}")
    N[3,3] = 1
    @test N[3,3] == one(BigFloat)

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
  print("""
Iteration end

""")

  print("""
Test similar
""")
  for ndim = 1:4
    @test size(similar(ResidueCount{Int, ndim, false}()).marginals) == (20, ndim)
    @test size(similar(ResidueCount{Int, ndim, true}()).marginals) == (21, ndim)
  end

  @test size(similar(ResidueCount{Int, 1, false}(),4).marginals) == (20, 4)
  @test isa(similar(ResidueCount{Int, 1, false}(),BigFloat).total, BigFloat)

  print("""
Test count! with Clusters
""")
  @test_throws BoundsError count!(zeros(ResidueCount{Float64, 1, true}), false_clusters, seq)
  @test count!(zeros(ResidueCount{Float64, 1, false}), false_clusters, seq).table == getweight(false_clusters)
  @test count!(zeros(ResidueCount{Float64, 2, false}), false_clusters, seq, seq).marginals[:,1] == getweight(false_clusters)
  @test count!(zeros(ResidueCount{Float64, 3, false}), false_clusters, seq, seq, seq).marginals[:,1] == getweight(false_clusters)

  print("""
Test count
""")
  @test count(seq, weight=false_clusters).table == getweight(false_clusters)
  @test count(seq, seq, weight=false_clusters).marginals[:,1] == getweight(false_clusters)
  @test count(seq, seq, seq, weight=false_clusters).marginals[:,1] == getweight(false_clusters)

  @test count(seq, usegap=true).total == 21
  @test count(seq, seq, usegap=true).total == 21
  @test count(seq, seq, seq, usegap=true).total == 21

  @test count(AdditiveSmoothing(1.0), seq, usegap=true).total == 21 + 21
  @test count(AdditiveSmoothing(1.0), seq, seq, usegap=true).total == 21 + (21*21)


  print("""

Tests for ResidueProbability
----------------------
""")

  @test size(ResidueProbability{BigFloat, 1, false}().table) == (20,)
  @test size(ResidueProbability{BigFloat, 1, true}().table) == (21,)

  @test size(ResidueProbability{BigFloat, 2, false}().table) == (20,20)
  @test size(ResidueProbability{BigFloat, 2, true}().table) == (21,21)

  @test size(ResidueProbability{BigFloat, 3, false}().table) == (20,20,20)
  @test size(ResidueProbability{BigFloat, 3, true}().table) == (21,21,21)

  for ndim = 1:4
    @test size(ResidueProbability{BigFloat, ndim, true}().marginals) == (21, ndim)
    @test size(ResidueProbability{BigFloat, ndim, false}().marginals) == (20, ndim)
  end

  print("""

Iterate for N & UseGap:
""")
  for ndim = 1:4, usegap = Bool[true, false]
    println("# Iteration = N: $(ndim) UseGap: $(usegap)")
    nres = usegap ? 21 : 20

    println("Test zeros for ResidueProbability{BigFloat, $(ndim), $(usegap)}")
    N = zeros(ResidueProbability{BigFloat, ndim, usegap})
    @test sum(N.table) == zero(BigFloat)

    println("Test iteration interface for ResidueProbability{BigFloat, $(ndim), $(usegap)}")
    @test collect(N) == zeros(BigFloat, nres^ndim)

    println("Test size for ResidueProbability{BigFloat, $(ndim), $(usegap)}")
    @test size(N) == size(N.table)

    println("Test indexing with i for ResidueProbability{BigFloat, $(ndim), $(usegap)}")
    @test N[1] == zero(BigFloat)
    @test N[Int(Residue('R'))] == zero(BigFloat)

    println("Test setindex! with i for ResidueProbability{BigFloat, $(ndim), $(usegap)}")
    N[3] = 1.0
    @test N[3] == one(BigFloat)

    println("Test update! for ResidueProbability{BigFloat, $(ndim), $(usegap)}")
    N = zeros(ResidueProbability{BigFloat, ndim, usegap})
    N[1] = 0.5
    update!(N)
    @test N.marginals[1,1] == big"0.5"

    println("Test normalize! for ResidueProbability{BigFloat, $(ndim), $(usegap)}")
    normalize!(N)
    @test N[1] == one(BigFloat)
    @test N.marginals[1,1] == one(BigFloat)
    @test sum(N) == one(BigFloat)

    println("Test fill! with ResidueCount{Int, $(ndim), $(usegap)} for ResidueProbability{BigFloat, $(ndim), $(usegap)}")
    C = fill!(zeros(ResidueCount{Int, ndim, usegap}), AdditiveSmoothing(1))
    P = fill!(N, C)
    @test P[1] == one(BigFloat) / length(C.table)
    @test P.marginals[1,1] == one(BigFloat) / nres
    @test_approx_eq_eps sum(P) one(BigFloat) 1e-72
  end

  print("""
Iteration end

""")

  print("""
Test Float64
""")

  let N = zeros(ResidueProbability{Float64, 3, true}), C = fill!(zeros(ResidueCount{Int, 3, true}), AdditiveSmoothing(1))
    P = fill!(N, C)
    @test_approx_eq sum(P) one(Float64)
    @test isa(sum(P), Float64)
  end

  print("""

Iterate for UseGap in [true, false]:
""")
  for usegap = Bool[true, false]
    println("# Iteration = UseGap: $(usegap)")
    nres = usegap ? 21 : 20

    N = zeros(ResidueProbability{BigFloat, 2, usegap})

    println("Test indexing with ij for ResidueProbability{BigFloat, 2, $(usegap)}")
    @test N[1,1] == zero(BigFloat)
    @test N[Int(Residue('R')), Int(Residue('R'))] == zero(BigFloat)

    println("Test setindex! with ij for ResidueProbability{BigFloat, 2, $(usegap)}")
    N[3,3] = 1.0
    @test N[3,3] == one(BigFloat)
  end
  print("""
Iteration end

""")

  for ndim = 1:4
    @test size(similar(ResidueProbability{BigFloat, ndim, false}()).marginals) == (20, ndim)
    @test size(similar(ResidueProbability{BigFloat, ndim, true}()).marginals) == (21, ndim)
  end

  @test size(similar(ResidueProbability{BigFloat, 1, false}(),4).marginals) == (20, 4)

  print("""
Test pseudofrequencies
""")
  Pab = fill!(zeros(ResidueProbability{BigFloat, 2, false}), count(seq, seq))

  Gab = blosum_pseudofrequencies!(zeros(ResidueProbability{BigFloat, 2, false}), Pab)
  @test_approx_eq sum(Gab) 1.0

  @test_approx_eq sum( apply_pseudofrequencies!(copy(Pab), Gab, 1, 1) ) 1.0

  Pab_with_pseudofrequency = apply_pseudofrequencies!(copy(Pab), Gab, 10, 0)
  @test Pab_with_pseudofrequency == Pab


  print("""
Test probabilities
""")

  @test probabilities(Float64, seq, weight=false_clusters).table == getweight(false_clusters)/sum(getweight(false_clusters))
  @test probabilities(Float64, seq, seq, weight=false_clusters).marginals[:,1] == getweight(false_clusters)/sum(getweight(false_clusters))
  @test probabilities(Float64, seq, seq, seq, weight=false_clusters).marginals[:,1] == getweight(false_clusters)/sum(getweight(false_clusters))


  @test_approx_eq probabilities(Float64, seq)[1] 0.05
  @test_approx_eq probabilities(Float32, seq, seq)[2,2] Float32(0.05)
  @test_approx_eq probabilities(BigFloat, seq, seq, seq)[3,3,3] one(BigFloat)/20

  @test_approx_eq probabilities(Float64, AdditiveSmoothing(1.0), seq, usegap=true)[1] 1/21.0
  Pab = probabilities(Float64, AdditiveSmoothing(1.0), seq, seq, usegap=true)
  @test_approx_eq Pab[1,1] 2.0/(21 + 21*21)
  @test_approx_eq Pab[1,2] 1.0/(21 + 21*21)
  @test_approx_eq sum(Pab) 1.0

  @test_approx_eq probabilities(Float64, 1, 0, seq, seq)[2,2] 0.05

  print("""
Test delete_dimensions
""")
  Pxyz = probabilities(Float64, seq, seq, seq);
  @test delete_dimensions(Pxyz, 3) == probabilities(Float64, seq, seq)
  @test delete_dimensions(Pxyz, 3, 2) == probabilities(Float64, seq)
  Pxy = delete_dimensions(Pxyz, 3);
  @test delete_dimensions!(Pxy, Pxyz, 1) == probabilities(Float64, seq, seq)
  @test_approx_eq sum(Pxy) 1.0
  Nxyz = count(seq, seq, seq);
  Nxy = delete_dimensions(Nxyz, 3);
  @test Nxy == count(seq, seq)
  @test delete_dimensions(Nxyz, 3, 2) == count(seq)
  @test delete_dimensions!(Nxy, Nxyz, 1) == count(seq, seq)
  @test sum(Nxy) == Nxyz.total

end
