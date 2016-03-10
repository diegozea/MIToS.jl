# using Base.Test
# using MIToS.MSA
# using MIToS.Information

let s = res"ARNDCQEGHILKMFPSTWYV-",
  r = reverse(s),
  g = res"GGGGGGGGGGGGGGGGGGGG",
  Pg = probabilities(BigFloat, g),
  Pgg = probabilities(BigFloat, g,g),
  Pggg = probabilities(BigFloat, g,g,g),
  Ps = probabilities(BigFloat, s),
  Pss = probabilities(BigFloat, s,s),
  Psr = probabilities(BigFloat, s,r),
  Psss = probabilities(BigFloat, s,s,s),
  prob_random = probabilities(BigFloat, AdditiveSmoothing(1.0), Residue[], Residue[]),
  count_random = count(AdditiveSmoothing(1.0), Residue[], Residue[])

  print("""

Tests for Information Measures
==============================
""")

  print("""

Entropy
-------
""")

  print("""
Entropy: H(X)
""")

  @test estimate(Entropy{BigFloat}(), Pg) == zero(BigFloat)
  @test_approx_eq estimate(Entropy{BigFloat}(), Ps) log(big"20")
  @test_approx_eq estimate(Entropy{BigFloat}(2), Ps) log(2, big"20")

  print("""
Joint Entropy: H(X,Y)
""")

  @test estimate(Entropy{BigFloat}(), Pgg) == zero(BigFloat)
  @test_approx_eq estimate(Entropy{BigFloat}(), Pss) log(big"20")
  @test_approx_eq estimate(Entropy{BigFloat}(), Psr) log(big"19")
  @test_approx_eq estimate(Entropy{BigFloat}(2), Pss) log(2, big"20")
  @test_approx_eq estimate(Entropy{BigFloat}(), count_random) log(big"400")
  @test_approx_eq estimate(Entropy{BigFloat}(),  prob_random) log(big"400")

  print("""
Joint Entropy: H(X,Y,Z)
""")

  @test estimate(Entropy{BigFloat}(), Pggg) == zero(BigFloat)
  @test_approx_eq estimate(Entropy{BigFloat}(), Psss) log(big"20")

  print("""
Entropy using ResidueCount
""")
  @test estimate(Entropy{BigFloat}(), count(g,g)) == estimate(Entropy{BigFloat}(), Pgg)
  @test estimate(Entropy{BigFloat}(), count(s,s)) == estimate(Entropy{BigFloat}(), Pss)
  @test_approx_eq estimate(Entropy{BigFloat}(), count(s,r)) estimate(Entropy{BigFloat}(), Psr)

  print("""

Mutual Information
------------------
""")

  @test estimate(MutualInformation{BigFloat}(), Pgg) == zero(BigFloat)
  @test estimate(MutualInformation{BigFloat}(), probabilities(BigFloat, g,s)) == zero(BigFloat)
  @test_approx_eq estimate(MutualInformation{BigFloat}(),  prob_random) zero(BigFloat) # diff for Float64 is 1e-15

  print("""
MI(X,Y) = H(X) + H(Y) - H(X,Y)
""")
  @test estimate(MutualInformation{BigFloat}(), Psr) == estimate_on_marginal(Entropy{BigFloat}(), Psr, 1) + estimate_on_marginal(Entropy{BigFloat}(), Psr, 2) - estimate(Entropy{BigFloat}(), Psr)
  @test estimate(MutualInformation{BigFloat}(), Pss) == estimate_on_marginal(Entropy{BigFloat}(), Pss, 1) + estimate_on_marginal(Entropy{BigFloat}(), Pss, 2) - estimate(Entropy{BigFloat}(), Pss)
  @test estimate(MutualInformation{BigFloat}(2), Psr) == estimate_on_marginal(Entropy{BigFloat}(2), Psr, 1) + estimate_on_marginal(Entropy{BigFloat}(2), Psr, 2) - estimate(Entropy{BigFloat}(2), Psr)
  @test estimate(MutualInformation{BigFloat}(20), Pss) == estimate_on_marginal(Entropy{BigFloat}(20), Pss, 1) + estimate_on_marginal(Entropy{BigFloat}(20), Pss, 2) - estimate(Entropy{BigFloat}(20), Pss)
  @test_approx_eq_eps estimate(MutualInformation{BigFloat}(), prob_random) ( estimate_on_marginal(Entropy{BigFloat}(), prob_random, 1) +
                                                                              estimate_on_marginal(Entropy{BigFloat}(),  prob_random, 2) -
                                                                              estimate(Entropy{BigFloat}(), prob_random) ) 1e-74 # diff for Float64 is 1e-15

  print("""
MI using ResidueCount
""")
  @test_approx_eq estimate(MutualInformation{BigFloat}(), count_random) zero(BigFloat)
  @test estimate(MutualInformation{BigFloat}(), count(g,g)) == estimate(MutualInformation{BigFloat}(), Pgg)
  @test estimate(MutualInformation{BigFloat}(), count(s,s)) == estimate(MutualInformation{BigFloat}(), Pss)
  @test_approx_eq estimate(MutualInformation{BigFloat}(), count(s,r)) estimate(MutualInformation{BigFloat}(), Psr)

  print("""
MI(X,Y,Z)
""")
  @test estimate(MutualInformation{BigFloat}(), Pggg) == zero(BigFloat)
  @test estimate(MutualInformation{BigFloat}(), probabilities(BigFloat, g,s,r)) == zero(BigFloat)
  @test_approx_eq estimate(MutualInformation{BigFloat}(), probabilities(BigFloat, s,s,r)) estimate(MutualInformation{BigFloat}(), count(s,s,r))

  print("""
MI(X,Y,Z) = H(X) + H(Y) + H(Z) - H(X,Y) - H(X,Z) - H(Y,Z) + H(X,Y,Z)
""")
  @test estimate(MutualInformation{BigFloat}(), Psss) == ( estimate_on_marginal(Entropy{BigFloat}(), Psss, 1) + estimate_on_marginal(Entropy{BigFloat}(), Psss, 2) + estimate_on_marginal(Entropy{BigFloat}(), Psss, 3) -
                                                            estimate(Entropy{BigFloat}(), Pss) - estimate(Entropy{BigFloat}(), Pss) - estimate(Entropy{BigFloat}(), Pss) + estimate(Entropy{BigFloat}(), Psss) )

  print("""
MI(X,Y,Z) <= min{ H(X,Y), H(X,Z), H(Y,Z) }
""")
  @test_approx_eq estimate(MutualInformation{BigFloat}(), Psss) estimate(MutualInformation{BigFloat}(),Pss)

  print("""

Pairwise Gap Percentage
-----------------------
""")

  @test estimate(GapUnionPercentage{Float64}(),        count(res"AA--", res"--AA", usegap=true)) == 100.0
  @test estimate(GapIntersectionPercentage{Float64}(), count(res"AA--", res"--AA", usegap=true)) == 0.0

  @test estimate(GapUnionPercentage{Float64}(),        count(res"AA--", res"--AA", usegap=true, weight=Float64[.25, .25, .25, .25])) == 100.0
  @test estimate(GapIntersectionPercentage{Float64}(), count(res"AA--", res"--AA", usegap=true, weight=Float64[.25, .25, .25, .25])) == 0.0

  @test estimate(GapUnionPercentage{Float64}(),        count(res"AAA-", res"AA--", usegap=true)) == 50.0
  @test estimate(GapIntersectionPercentage{Float64}(), count(res"AAA-", res"AA--", usegap=true)) == 25.0

  @test estimate(GapUnionPercentage{Float64}(),        count(res"AAA-", res"AA--", usegap=true, weight=Float64[.2, .2, .2, .4])) == 60.0
  @test estimate(GapIntersectionPercentage{Float64}(), count(res"AAA-", res"AA--", usegap=true, weight=Float64[.2, .2, .2, .4])) == 40.0
end
