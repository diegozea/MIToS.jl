abstract InformationMeasure

abstract SymmetricMeasure <: InformationMeasure

#estimate_on_marginals(measure::InformationMeasure, table::ResidueContingencyTables, marginal::Int) = estimate(measure, table.marginal[:,marginal])
#estimate_on_marginals(measure::InformationMeasure, table::ResidueContingencyTables, marginal::Int, base::Real) = estimate(measure, table.marginal[:,marginal], base)

# Entropy

type Entropy <: SymmetricMeasure end

## Estimate Entropy using ResidueProbability

"""```estimate(Entropy(), p [, base])```

`p` should be a `ResidueProbability` table.
"""
function estimate{T, N, UseGap}(measure::Entropy, p::ResidueProbability{T, N, UseGap})
  H = zero(T)
  for i in 1:length(p)
    @inbounds pi = p[i]
    if pi != 0.0
      H += pi * log(pi)
    end
  end
  -H
end

"""```estimate_on_marginal(Entropy(), p, marginal[, base])```

This function estimate the entropy H(X) if marginal is 1, H(Y) for 2, etc.
"""
function estimate_on_marginal{T, N, UseGap}(measure::Entropy, p::ResidueProbability{T, N, UseGap}, marginal::Int)
  H = zero(T)
  for i in 1:nresidues(p)
    @inbounds pi = p.marginals[i, marginal]
    if pi != 0.0
      H += pi * log(pi)
    end
  end
  -H
end

## Estimate Entropy using ResidueCount

"""```estimate(Entropy(), n::ResidueCount [, base])```

It's the fastest option (you don't spend time on probability calculations)."""
function estimate(measure::Entropy, n::ResidueCount)
  H = zero(Float64)
  total = Float64(n.total)
  for i in 1:length(n)
    @inbounds ni = n[i]
    if ni != 0.0
      H += ni * log(ni/total)
    end
  end
  -H/total
end

function estimate_on_marginal(measure::Entropy, n::ResidueCount, marginal::Int)
  H = zero(Float64)
  total = Float64(n.total)
  for i in 1:nresidues(n)
    @inbounds ni = n.marginals[i, marginal]
    if ni != 0.0
      H += ni * log(ni/total)
    end
  end
  -H/total
end


estimate(measure::Entropy, p::ResidueContingencyTables, base::Real) = estimate(measure, p) / log(base)
estimate_on_marginal(measure::Entropy, table::ResidueContingencyTables, marginal::Int, base::Real) = estimate_on_marginal(measure, table, marginal) / log(base)

# Mutual Information

type MutualInformation <: SymmetricMeasure end


"""```estimate(MutualInformation(), pxy::ResidueProbability [, base])```

Calculate Mutual Information from `ResidueProbability`."""
function estimate{T, UseGap}(measure::MutualInformation, pxy::ResidueProbability{T, 2,UseGap})
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
  MI
end

estimate(measure::MutualInformation, pxy, base::Real) = estimate(measure, pxy) / log(base)

"""```estimate(MutualInformation(), pxy::ResidueCount [, base])```

Calculate Mutual Information from `ResidueCount`.
It's the fastest option (you don't spend time on probability calculations)."""
function estimate{T, UseGap}(measure::MutualInformation, nxy::ResidueCount{T, 2,UseGap})
  MI = zero(Float64)
  N = Float64(nxy.total)
  @inbounds for j in 1:nresidues(nxy)
    nj = nxy.marginals[j,2]
    if nj > 0.0
      for i in 1:nresidues(nxy)
        ni = nxy.marginals[i,1]
        nij = nxy[i,j]
        if nij > 0.0 && ni > 0.0
          MI +=  nij * log((N * nij)/(ni*nj))
        end
      end
    end
  end
  MI/N
end

estimate(measure::MutualInformation, pxy, base::Real) = estimate(measure, pxy) / log(base)

function estimate{T, UseGap}(measure::MutualInformation, pxyz::ResidueContingencyTables{T, 3, UseGap})
  pxy = delete_dimensions(pxyz,3)
  return( estimate_on_marginal(Entropy(), pxyz, 1) + # H(X)
  estimate_on_marginal(Entropy(), pxyz ,2) + # H(Y)
  estimate_on_marginal(Entropy(), pxyz, 3) - # H(Z)
  estimate(Entropy(), pxy) - # H(X, Y)
  estimate(Entropy(), delete_dimensions!(pxy, pxyz, 2)) - # H(X, Z)
  estimate(Entropy(), delete_dimensions!(pxy, pxyz, 1)) + # H(Y, Z)
  estimate(Entropy(), pxyz) ) # H(X, Y, Z)
end
