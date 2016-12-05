# Pseudocounts
# ============

abstract Pseudocount{T<:Real}

immutable NoPseudocount <: Pseudocount end

"""
**Additive Smoothing** or fixed pseudocount  `λ`  for `ResidueCount`
(in order to estimate probabilities when the number of samples is low).

Common values of `λ` are:

- `0` :  No cell frequency prior, gives you the maximum likelihood estimator.
- `0.05` is the optimum value for `λ` found in Buslje et. al. 2009, similar results was obtained for `λ` in the range [0.025, 0.075].
- `1 / p` : Perks prior (Perks, 1947) where `p` the number of parameters (i.e. residues, pairs of residues) to estimate. If `p` is the number of residues (`20` without counting gaps), this gives you `0.05`.
- `sqrt(n) / p` : Minimax prior (Trybula, 1958) where `n` is the number of samples and `p` the number of parameters to estimate.  If the number of samples `n` is 400 (minimum number of sequence clusters for achieve good performance in Buslje et. al. 2009) for estimating 400 parameters (pairs of residues without counting gaps) this gives you `0.05`.
- `0.5` : Jeffreys prior (Jeffreys, 1946).
- `1` : Bayes-Laplace uniform prior, aka. Laplace smoothing.
"""
immutable AdditiveSmoothing{T} <: Pseudocount{T}
  λ::T
end

Base.zero{T}(::Type{AdditiveSmoothing{T}}) = AdditiveSmoothing(zero(T))
Base.one{T}(::Type{AdditiveSmoothing{T}}) = AdditiveSmoothing(one(T))
