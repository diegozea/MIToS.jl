using Base.Test

using MIToS.MSA
const s = Residue(MSA._to_char);
const r = reverse(s);
const g = Residue("GGGGGGGGGGGGGGGGGGGG".data);

using MIToS.Information
const Pg = probabilities(g)
const Pgg = probabilities(g,g)
const Pggg = probabilities(g,g,g)
const Ps = probabilities(s)
const Pss = probabilities(s,s)
const Psr = probabilities(s,r)
const Psss = probabilities(s,s,s)

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

@test estimate(Entropy(), Pg) == 0.0
@test estimate(Entropy(), Ps) == log(20)
@test estimate(Entropy(), Ps, 2) == log(2, 20)

print("""
Joint Entropy: H(X,Y)
""")

@test estimate(Entropy(), Pgg) == 0.0
@test estimate(Entropy(), Pss) == log(20)
@test estimate(Entropy(), Psr) == log(19)
@test estimate(Entropy(), Pss, 2) == log(2, 20)

print("""
Joint Entropy: H(X,Y,Z)
""")

@test estimate(Entropy(), Pggg) == 0.0
@test estimate(Entropy(), Psss) == log(20)

print("""
Entropy using ResidueCount
""")
@test estimate(Entropy(), count(g,g)) == estimate(Entropy(), Pgg)
@test estimate(Entropy(), count(s,s)) == estimate(Entropy(), Pss)
@test_approx_eq estimate(Entropy(), count(s,r)) estimate(Entropy(), Psr)

print("""

Mutual Information
------------------
""")

@test estimate(MutualInformation(), Pgg) == 0.0
@test estimate(MutualInformation(), probabilities(g,s)) == 0.0

print("""
MI(X,Y) = H(X) + H(Y) - H(X,Y)
""")
@test estimate(MutualInformation(), Psr) == estimate_on_marginal(Entropy(), Psr, 1) + estimate_on_marginal(Entropy(), Psr, 2) - estimate(Entropy(), Psr)
@test estimate(MutualInformation(), Pss) == estimate_on_marginal(Entropy(), Pss, 1) + estimate_on_marginal(Entropy(), Pss, 2) - estimate(Entropy(), Pss)
@test estimate(MutualInformation(), Psr, 2) == estimate_on_marginal(Entropy(), Psr, 1, 2) + estimate_on_marginal(Entropy(), Psr, 2, 2) - estimate(Entropy(), Psr, 2)
@test estimate(MutualInformation(), Pss, 20) == estimate_on_marginal(Entropy(), Pss, 1, 20) + estimate_on_marginal(Entropy(), Pss, 2, 20) - estimate(Entropy(), Pss, 20)

print("""
MI using ResidueCount
""")
@test estimate(MutualInformation(), count(g,g)) == estimate(MutualInformation(), Pgg)
@test estimate(MutualInformation(), count(s,s)) == estimate(MutualInformation(), Pss)
@test_approx_eq estimate(MutualInformation(), count(s,r)) estimate(MutualInformation(), Psr)

print("""
MI(X,Y,Z)
""")
@test estimate(MutualInformation(), Pggg) == 0.0
@test estimate(MutualInformation(), probabilities(g,s,r)) == 0.0
@test_approx_eq_eps estimate(MutualInformation(), probabilities(AdditiveSmoothing(1e10),s,s,s)) 0.0 1e-10

print("""
MI(X,Y,Z) = H(X) + H(Y) + H(Z) - H(X,Y) - H(X,Z) - H(Y,Z) + H(X,Y,Z)
""")
@test estimate(MutualInformation(), Psss) == ( estimate_on_marginal(Entropy(), Psss, 1) + estimate_on_marginal(Entropy(), Psss, 2) + estimate_on_marginal(Entropy(), Psss, 3) -
  estimate(Entropy(), Pss) - estimate(Entropy(), Pss) - estimate(Entropy(), Pss) + estimate(Entropy(), Psss) )

print("""
MI(X,Y,Z) <= min{ H(X,Y), H(X,Z), H(Y,Z) }
""")
@test_approx_eq estimate(MutualInformation(), Psss) estimate(MutualInformation(),Pss)
