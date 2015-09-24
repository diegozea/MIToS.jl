abstract InformationMeasure{T}

abstract SymmetricMeasure{T} <: InformationMeasure{T}

#estimate_on_marginals(measure::InformationMeasure, table::ResidueContingencyTables, marginal::Int) = estimate(measure, table.marginal[:,marginal])
#estimate_on_marginals(measure::InformationMeasure, table::ResidueContingencyTables, marginal::Int, base::Real) = estimate(measure, table.marginal[:,marginal], base)

# Entropy
# =======

immutable Entropy{T} <: SymmetricMeasure{T}
  base::T
end

call{T}(::Type{Entropy{T}}) = Entropy(T(Base.e))

## Estimate Entropy using ResidueProbability

"""```estimate(Entropy(base), p)```

`p` should be a `ResidueProbability` table. The result type is determined by `base`.
"""
function estimate{B, T, N, UseGap}(measure::Entropy{B}, p::ResidueProbability{T, N, UseGap})
  H = zero(B)
  for i in 1:length(p)
    @inbounds pi = p[i]
    if pi != 0.0
      H += pi * log(pi)
    end
  end
  -H/log(measure.base)
end

"""```estimate_on_marginal(Entropy(base), p, marginal)```

This function estimate the entropy H(X) if marginal is 1, H(Y) for 2, etc.
The result type is determined by `base`.
"""
function estimate_on_marginal{B, T, N, UseGap}(measure::Entropy{B}, p::ResidueProbability{T, N, UseGap}, marginal::Int)
  H = zero(B)
  for i in 1:nresidues(p)
    @inbounds pi = p.marginals[i, marginal]
    if pi != 0.0
      H += pi * log(pi)
    end
  end
  -H/log(measure.base)
end

## Estimate Entropy using ResidueCount

"""```estimate(Entropy(), n::ResidueCount [, base])```

It's the fastest option (you don't spend time on probability calculations).
The result type is determined by the `base`."""
function estimate{T}(measure::Entropy{T}, n::ResidueCount)
  H = zero(T)
  total = T(n.total)
  for i in 1:length(n)
    @inbounds ni = n[i]
    if ni != 0.0
      H += ni * log(ni/total)
    end
  end
  (-H/total)/log(measure.base)
end

function estimate_on_marginal{T}(measure::Entropy{T}, n::ResidueCount, marginal::Int)
  H = zero(T)
  total = T(n.total)
  for i in 1:nresidues(n)
    @inbounds ni = n.marginals[i, marginal]
    if ni != 0.0
      H += ni * log(ni/total)
    end
  end
  (-H/total)/log(measure.base)
end

# Mutual Information
# ==================

immutable MutualInformation{T} <: SymmetricMeasure{T}
  base::T
end

call{T}(::Type{MutualInformation{T}}) = MutualInformation(T(Base.e))

"""```estimate(MutualInformation(), pxy::ResidueProbability [, base])```

Calculate Mutual Information from `ResidueProbability`. The result type is determined by `base`."""
function estimate{B, T, UseGap}(measure::MutualInformation{B}, pxy::ResidueProbability{T, 2,UseGap})
  MI = zero(T)
  @inbounds for j in 1:nresidues(pxy)
    pj = pxy.marginals[j,2]
    if pj > 0.0
      for i in 1:nresidues(pxy)
        pi = pxy.marginals[i,1]
        pij = pxy[i,j]
        if pij > 0.0 && pi > 0.0
          MI +=  pij * log(pij/(pi*pj))
        end
      end
    end
  end
  MI/log(measure.base)
end

@inline _mi{T}(N::T, nij, ni, nj) = ifelse(nij > zero(T) && ni > zero(T), T(nij * log((N * nij)/(ni * nj))), zero(T))

"""```estimate(MutualInformation(), pxy::ResidueCount [, base])```

Calculate Mutual Information from `ResidueCount`. The result type is determined by the `base`.
It's the fastest option (you don't spend time on probability calculations)."""
function estimate{B, T, UseGap}(measure::MutualInformation{B}, nxy::ResidueCount{T, 2,UseGap})
  MI = zero(B)
  N = B(nxy.total)
  marginals = nxy.marginals
  @inbounds for j in 1:nresidues(nxy)
    nj = marginals[j,2]
    if nj > 0.0
      @inbounds @simd for i in 1:nresidues(nxy)
        MI += _mi(N, nxy[i,j], marginals[i,1], nj)
#         ni = nxy.marginals[i,1]
#         nij = nxy[i,j]
#         if nij > 0.0 && ni > 0.0
#           MI +=  nij * log((N * nij)/(ni*nj))
#         end
      end
    end
  end
  (MI/N)/log(measure.base)
end

function estimate{B, T, UseGap}(measure::MutualInformation{B}, pxyz::ResidueContingencyTables{T, 3, UseGap})
  pxy = delete_dimensions(pxyz,3)
  return( estimate_on_marginal(Entropy(measure.base), pxyz, 1) + # H(X)
  estimate_on_marginal(Entropy(measure.base), pxyz ,2) + # H(Y)
  estimate_on_marginal(Entropy(measure.base), pxyz, 3) - # H(Z)
  estimate(Entropy(measure.base), pxy) - # H(X, Y)
  estimate(Entropy(measure.base), delete_dimensions!(pxy, pxyz, 2)) - # H(X, Z)
  estimate(Entropy(measure.base), delete_dimensions!(pxy, pxyz, 1)) + # H(Y, Z)
  estimate(Entropy(measure.base), pxyz) ) # H(X, Y, Z)
end

# Normalized Mutual Information by Entropy
# ========================================

# nMI(X, Y) = MI(X, Y) / H(X, Y)

immutable MutualInformationOverEntropy{T} <: SymmetricMeasure{T}
  base::T
end

call{T}(::Type{MutualInformationOverEntropy{T}}) = MutualInformationOverEntropy(T(Base.e))

function estimate{B}(measure::MutualInformationOverEntropy{B}, table)
  H = estimate(Entropy(measure.base), table)
  if H != zero(B)
    MI = estimate(MutualInformation(measure.base), table)
    return(MI/H)
  else
    return(zero(B))
  end
end
