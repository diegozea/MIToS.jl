using Base.Test

using MIToS.MSA
const s = res"ARNDCQEGHILKMFPSTWYV-"
const r = reverse(s);
const g = res"GGGGGGGGGGGGGGGGGGGG";

using MIToS.Information
const Pg = probabilities(BigFloat, g)
const Pgg = probabilities(BigFloat, g,g)
const Pggg = probabilities(BigFloat, g,g,g)
const Ps = probabilities(BigFloat, s)
const Pss = probabilities(BigFloat, s,s)
const Psr = probabilities(BigFloat, s,r)
const Psss = probabilities(BigFloat, s,s,s)

const  prob_random = probabilities(BigFloat, AdditiveSmoothing(1.0), Residue[], Residue[])
const count_random = count(AdditiveSmoothing(1.0), Residue[], Residue[])

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

@test estimate(Entropy(), Pg) == zero(BigFloat)
@test_approx_eq estimate(Entropy(), Ps) log(big"20")
@test_approx_eq_eps estimate(Entropy(), Ps, 2) log(2, big"20") 1e-15

print("""
Joint Entropy: H(X,Y)
""")

@test estimate(Entropy(), Pgg) == zero(BigFloat)
@test_approx_eq estimate(Entropy(), Pss) log(big"20")
@test_approx_eq estimate(Entropy(), Psr) log(big"19")
@test_approx_eq_eps estimate(Entropy(), Pss, 2) log(2, big"20") 1e-15
@test_approx_eq estimate(Entropy(), count_random) log(big"400")
@test_approx_eq estimate(Entropy(),  prob_random) log(big"400")

print("""
Joint Entropy: H(X,Y,Z)
""")

@test estimate(Entropy(), Pggg) == zero(BigFloat)
@test_approx_eq estimate(Entropy(), Psss) log(big"20")

print("""
Entropy using ResidueCount
""")
@test estimate(Entropy(), count(g,g)) == estimate(Entropy(), Pgg) # Entropy & count is Float64, but zero(Float64) == zero(BigFloat)
@test estimate(Entropy(), count(s,s)) == Float64(estimate(Entropy(), Pss)) # Entropy & count is Float64
@test_approx_eq estimate(Entropy(), count(s,r)) estimate(Entropy(), Psr)

print("""

Mutual Information
------------------
""")

@test estimate(MutualInformation(), Pgg) == zero(BigFloat)
@test estimate(MutualInformation(), probabilities(BigFloat, g,s)) == zero(BigFloat)
@test_approx_eq estimate(MutualInformation(),  prob_random) zero(BigFloat) # eps for Float64 is 1e-15

print("""
MI(X,Y) = H(X) + H(Y) - H(X,Y)
""")
@test estimate(MutualInformation(), Psr) == estimate_on_marginal(Entropy(), Psr, 1) + estimate_on_marginal(Entropy(), Psr, 2) - estimate(Entropy(), Psr)
@test estimate(MutualInformation(), Pss) == estimate_on_marginal(Entropy(), Pss, 1) + estimate_on_marginal(Entropy(), Pss, 2) - estimate(Entropy(), Pss)
@test estimate(MutualInformation(), Psr, 2) == estimate_on_marginal(Entropy(), Psr, 1, 2) + estimate_on_marginal(Entropy(), Psr, 2, 2) - estimate(Entropy(), Psr, 2)
@test estimate(MutualInformation(), Pss, 20) == estimate_on_marginal(Entropy(), Pss, 1, 20) + estimate_on_marginal(Entropy(), Pss, 2, 20) - estimate(Entropy(), Pss, 20)
@test_approx_eq_eps estimate(MutualInformation(), prob_random) ( estimate_on_marginal(Entropy(), prob_random, 1) +
                                                                  estimate_on_marginal(Entropy(),  prob_random, 2) -
                                                                  estimate(Entropy(), prob_random) ) 1e-74 # eps for Float64 is 1e-15

print("""
MI using ResidueCount
""")
@test_approx_eq estimate(MutualInformation(), count_random) zero(BigFloat)
@test estimate(MutualInformation(), count(g,g)) == estimate(MutualInformation(), Pgg) # MI & count is Float64, but zero(Float64) == zero(BigFloat)
@test estimate(MutualInformation(), count(s,s)) == Float64(estimate(MutualInformation(), Pss)) # MI & count is Float64
@test_approx_eq estimate(MutualInformation(), count(s,r)) estimate(MutualInformation(), Psr)

print("""
MI(X,Y,Z)
""")
@test estimate(MutualInformation(), Pggg) == zero(BigFloat)
@test estimate(MutualInformation(), probabilities(BigFloat, g,s,r)) == zero(BigFloat)
@test_approx_eq_eps estimate(MutualInformation(), probabilities(BigFloat, AdditiveSmoothing(1e10),s,s,s)) zero(BigFloat) 1e-22 # eps for Float64 1e-10
@test_approx_eq estimate(MutualInformation(), probabilities(BigFloat, s,s,r)) estimate(MutualInformation(), count(s,s,r))

print("""
MI(X,Y,Z) = H(X) + H(Y) + H(Z) - H(X,Y) - H(X,Z) - H(Y,Z) + H(X,Y,Z)
""")
@test estimate(MutualInformation(), Psss) == ( estimate_on_marginal(Entropy(), Psss, 1) + estimate_on_marginal(Entropy(), Psss, 2) + estimate_on_marginal(Entropy(), Psss, 3) -
  estimate(Entropy(), Pss) - estimate(Entropy(), Pss) - estimate(Entropy(), Pss) + estimate(Entropy(), Psss) )

print("""
MI(X,Y,Z) <= min{ H(X,Y), H(X,Z), H(Y,Z) }
""")
@test_approx_eq estimate(MutualInformation(), Psss) estimate(MutualInformation(),Pss)
