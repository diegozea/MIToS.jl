import Base: zero, zeros, copy, deepcopy, fill!, getindex, setindex!

### PROBABILITIES AND PSEUDOCOUNTS ###

abstract Pseudocount

immutable Fixed <: Pseudocount
  λ::Float64
end

type Pseudofrequencies <: Pseudocount
  α::Float64
  β::Float64
  Gab::Matrix{Float64}
end

type ResidueProbabilities
  Pa::Vector{Float64}
end

type ResiduePairProbabilities
  Pab::Matrix{Float64}
  Pa::Vector{Float64}
  Pb::Vector{Float64}
end

### Constructors ###

zero(::Type{Fixed}) = Fixed(zero(Float64))

zeros(::Type{Pseudofrequencies}) = Pseudofrequencies(zero(Float64), zero(Float64), zeros(Float64, (20,20)))
zeros(::Type{ResidueProbabilities}) = ResidueProbabilities(zeros(Float64,20))
zeros(::Type{ResiduePairProbabilities}) = ResiduePairProbabilities(zeros(Float64, (20,20)), zeros(Float64,20), zeros(Float64,20))

Fixed() = zero(Fixed)
Pseudofrequencies() = Pseudofrequencies(zero(Float64), zero(Float64), Array(Float64, (20,20)))
ResidueProbabilities() = ResidueProbabilities(Array(Float64,20))
ResiduePairProbabilities() = ResiduePairProbabilities(Array(Float64, (20,20)), Array(Float64,20), Array(Float64,20))

Pseudofrequencies(α::Float64, β::Float64) = Pseudofrequencies(α, β, Array(Float64, (20,20)))

### Copy ###

copy(x::Fixed) = Fixed(copy(x.λ))
copy(x::Pseudofrequencies) = Pseudofrequencies(copy(x.α), copy(x.β), copy(x.Gab))
copy(pa::ResidueProbabilities) = ResidueProbabilities(copy(pa.Pa))
copy(pab::ResiduePairProbabilities) = ResiduePairProbabilities(copy(pab.Pab), copy(pab.Pa), copy(pab.Pb))

### Indexing ###

getindex(pa::ResidueProbabilities, i::Int) = getindex(pa.Pa, i)
setindex!(pa::ResidueProbabilities, p::Float64, i::Int) = setindex!(pa.Pa, p, i)

getindex(pab::ResiduePairProbabilities, i::Int, j::Int) = getindex(pab.Pab, i, j)

getindex(pse::Pseudofrequencies, i::Int, j::Int) = getindex(pse.Gab, i, j)

### Estimation (Filling the Probabilities) ###

"""Fill the Pab matrix with λ, but doesn't update the marginal probabilities"""
function __initialize!(pab::ResiduePairProbabilities)
  fill!(pab.Pab, zero(Float64))
  zero(Float64)
end

function __initialize!(pab::ResiduePairProbabilities, pseudocount::Fixed)
  fill!(pab.Pab, pseudocount.λ)
  pseudocount.λ * 400.0
end

function __finalize!(pab::ResiduePairProbabilities, total::Float64)
  pab.Pab[:,:] /= total
  pab.Pa[:] = sum(pab.Pab,1)
  pab.Pa[:] = sum(pab.Pab,2)
  pab
end

function __fill_kernel!(pab::ResiduePairProbabilities, seqi::AbstractVector{Residue}, seqj::AbstractVector{Residue}, total::Float64)
  for i in 1:length(seqi) # seqj should have the same length
    if seqi[i] != GAP  && seqj[i] != GAP
      @inbounds pab.Pab[ seqi[i] , seqj[i] ] += one(Float64)
      total += one(Float64)
    end
  end
  total
end

function __fill_kernel!(pab::ResiduePairProbabilities, cl::Clusters, seqi::AbstractVector{Residue}, seqj::AbstractVector{Residue}, total::Float64)
  for i in 1:length(seqi) # seqj should have the same length
    if seqi[i] != GAP  && seqj[i] != GAP
      weight = getweight(cl, i)
      @inbounds pab.Pab[ seqi[i] , seqj[i] ] += weight
      total += weight
    end
  end
  total
end

function __fill_pseudofrequencies!(pseudocount::Pseudofrequencies, pab::ResiduePairProbabilities)
  total = zero(Float64)
  for a in 1:20, b in 1:20
    @inbounds pseudocount.Gab[a,b] = zero(Float64)
    for i in 1:20, j in 1:20
      P = pab[i,j]
      if P != 0
	      # BLOSUM62_P_i_j[i,a] is p(a | i)
	      @inbounds val = P * BLOSUM62_Pij[i,a] * BLOSUM62_Pij[j,b]
	      @inbounds pseudocount.Gab[a,b] += val
	      total += val
      end
    end
  end
  pseudocount.Gab[:,:] /= total
  pseudocount
end

function __apply_pseudofrequencies!(pab::ResiduePairProbabilities, pseudocount::Pseudofrequencies)
  frac = one(Float64) / ( pseudocount.α + pseudocount.β )
  total = zero(Float64)
  for i in 1:20, j in 1:20
    @inbounds val = ( pseudocount.α * pab[i,j] + pseudocount.β * pseudocount[i,j] ) * frac
    @inbounds pab.Pab[i,j] = val
    total += val
  end
  total
end

function fill!(pab::ResiduePairProbabilities, seqi::AbstractVector{Residue}, seqj::AbstractVector{Residue})
  total = __initialize!(pab)
  total = __fill_kernel!(pab, seqi, seqj, total)
  __finalize!(pab, total)
end

function fill!(pab::ResiduePairProbabilities, pseudocount::Fixed, seqi::AbstractVector{Residue}, seqj::AbstractVector{Residue})
  total = __initialize!(pab, pseudocount)
  total = __fill_kernel!(pab, seqi, seqj, total)
  __finalize!(pab, total)
end

function fill!(pab::ResiduePairProbabilities, pseudocount::Fixed, cl::Clusters, seqi::AbstractVector{Residue}, seqj::AbstractVector{Residue})
  total = __initialize!(pab, pseudocount)
  total = __fill_kernel!(pab, cl, seqi, seqj, total)
  __finalize!(pab, total)
end

function fill!(pab::ResiduePairProbabilities, pseudocount::Pseudofrequencies, seqi::AbstractVector{Residue}, seqj::AbstractVector{Residue})
  total = __initialize!(pab)
  total = __fill_kernel!(pab, seqi, seqj, total)
  if pseudocount.β != zero(Float64)
    __fill_pseudofrequencies!(pseudocount, pab)
    total = __apply_pseudofrequencies(pab, pseudocount)
  end
  __finalize!(pab, total)
end

function fill!(pab::ResiduePairProbabilities, pseudocount::Pseudofrequencies, cl::Clusters, seqi::AbstractVector{Residue}, seqj::AbstractVector{Residue})
  total = __initialize!(pab)
  total = __fill_kernel!(pab, cl, seqi, seqj, total)
  if pseudocount.β != zero(Float64)
    __fill_pseudofrequencies!(pseudocount, pab)
    total = __apply_pseudofrequencies(pab, pseudocount)
  end
  __finalize!(pab, total)
end
