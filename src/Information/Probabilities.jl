import  Base: zero, one, zeros, start, next, done, length, eltype,
        size, setindex!, getindex, similar, fill!, count #, print # , copy, deepcopy, fill!, getindex, setindex!

# Counts and Pseudocounts

abstract Pseudocount

"""
**Additive Smoothing** or fixed pseudocount  `λ`  for `ResidueCount` (in order to estimate probabilities when the number of samples is low).

Common values of `λ` are:

- `0` :  No cell frequency prior, gives you the maximum likelihood estimator.
- `0.05` is the optimum value for `λ` found in Buslje et. al. 2009, similar results was obtained for `λ` in the range [0.025, 0.075].
- `1 / p` : Perks prior (Perks, 1947) where `p` the number of parameters (i.e. residues, pairs of residues) to estimate. If `p` is the number of residues (`20` without counting gaps), this gives you `0.05`.
- `sqrt(n) / p` : Minimax prior (Trybula, 1958) where `n` is the number of samples and `p` the number of parameters to estimate.  If the number of samples `n` is 400 (minimum number of sequence clusters for achieve good performance in Buslje et. al. 2009) for estimating 400 parameters (pairs of residues without counting gaps) this gives you `0.05`.
- `0.5` : Jeffreys prior (Jeffreys, 1946).
- `1` : Bayes-Laplace uniform prior, aka. Laplace smoothing."""
immutable AdditiveSmoothing{T<:Real} <: Pseudocount
  λ::T
end

zero{T}(::Type{AdditiveSmoothing{T}}) = AdditiveSmoothing{T}(zero(T))
one{T}(::Type{AdditiveSmoothing{T}}) = AdditiveSmoothing{T}(one(T))

## Counts

"""`ResidueCount{T, N, UseGap}` is used for counting residues in columns (or sequences) of an MSA.
`N` is the dimensionality and should be an `Int`, i.e. 2 if 2 columns are used for counting pairs.
`UseGap` is a `Bool`, `true` means that **ResidueCount** counts gaps in the position 21."""
# The field marginal is used for pre allocation of marginal sums.
# The field total is used for saving the total.
type ResidueCount{T, N, UseGap} <: AbstractArray{T, N}
  counts::Array{T, N}
  marginals::Array{T, 2}
  total::T
end

call{T}(::Type{ResidueCount{T, 1, true}}) = ResidueCount{T, 1, true}(Array(T, 21), Array(T, (21,1)), zero(T))
call{T}(::Type{ResidueCount{T, 2, true}}) = ResidueCount{T, 2, true}(Array(T, (21, 21)), Array(T, (21,2)), zero(T))

call{T}(::Type{ResidueCount{T, 1, false}}) = ResidueCount{T, 1, false}(Array(T, 20), Array(T, (20,1)), zero(T))
call{T}(::Type{ResidueCount{T, 2, false}}) = ResidueCount{T, 2, false}(Array(T, (20, 20)), Array(T, (20,2)), zero(T))

function call{T, N, UseGap}(::Type{ResidueCount{T, N, UseGap}})
  nres = UseGap ? 21 : 20
  ResidueCount{T, N, UseGap}(Array(T, (Int[ nres for d in 1:N]...)), Array(T, (nres, N)), zero(T))
end

zeros{T}(::Type{ResidueCount{T, 1, true}}) = ResidueCount{T, 1, true}(zeros(T, 21), zeros(T, (21, 1)), zero(T))
zeros{T}(::Type{ResidueCount{T, 2, true}}) = ResidueCount{T, 2, true}(zeros(T, (21, 21)), zeros(T, (21, 2)), zero(T))

zeros{T}(::Type{ResidueCount{T, 1, false}}) = ResidueCount{T, 1, false}(zeros(T, 20), zeros(T, (20, 1)), zero(T))
zeros{T}(::Type{ResidueCount{T, 2, false}}) = ResidueCount{T, 2, false}(zeros(T, (20, 20)), zeros(T, (20, 2)), zero(T))

function zeros{T, N, UseGap}(::Type{ResidueCount{T, N, UseGap}})
  nres = UseGap ? 21 : 20
  ResidueCount{T, N, UseGap}(zeros(T, (Int[ nres for d in 1:N]...)), zeros(T, (nres, N)), zero(T))
end

### Abstract Array Interface

Base.linearindexing(::ResidueCount)	= Base.LinearFast()

"""`getindex(n::ResidueCount, i::Int)` gives you access to the counts, use `Int()` for indexing with `Residues`"""
getindex(n::ResidueCount, i::Int)	= getindex(n.counts, i)

"""`setindex!(n::ResidueCount, v, i::Int)` set a value into counts field, but doesn't update the marginals and total.
Use `Int()` for indexing with `Residues`. Use `update!` in order to calculate again the marginals and total."""
setindex!(n::ResidueCount, v, i::Int) = setindex!(n.counts, v, i)

#### Length & Size

length{T}(n::ResidueCount{T, 1, true}) = 21
length{T}(n::ResidueCount{T, 1,false}) = 20

length{T}(n::ResidueCount{T, 2, true}) = 441
length{T}(n::ResidueCount{T, 2,false}) = 400

length{T,N,UseGap}(n::ResidueCount{T, N, UseGap}) = length(n.counts)

size{T}(n::ResidueCount{T, 1, true}) = (21,)
size{T}(n::ResidueCount{T, 1,false}) = (20,)

size{T}(n::ResidueCount{T, 2, true}) = (21, 21)
size{T}(n::ResidueCount{T, 2,false}) = (20, 20)

size{T,N,UseGap}(n::ResidueCount{T, N, UseGap}) = size(n.counts)

#### Iteration Interface

start(n::ResidueCount) = 1

@inbounds next(n::ResidueCount, state::Int) = (n.counts[state], state + 1)

done{T,N,UseGap}(n::ResidueCount{T, N, UseGap}, state) = state > length(n)

eltype{T,N,UseGap}(::Type{ResidueCount{T, N, UseGap}}) = T

#### Similar

similar{T,N,UseGap}(n::ResidueCount{T,N,UseGap}) = ResidueCount{T,N,UseGap}()
similar{T,N,S,UseGap}(n::ResidueCount{T,N,UseGap}, ::Type{S}) = ResidueCount{S,N,UseGap}()
similar{T,N,UseGap}(n::ResidueCount{T,N,UseGap}, D::Int) = ResidueCount{T,D,UseGap}()
similar{T,N,S,UseGap}(n::ResidueCount{T,N,UseGap}, ::Type{S}, D::Int) = ResidueCount{S,D,UseGap}()

### Update!

function _tuple_without_index(index::Int, N::Int)
	used = Array(Int, N - 1)
	j = 1
	for i in 1:N
		if i != index
			used[j] = i
			j += 1
		end
	end
	(used...)
end

function update!{T,UseGap}(n::ResidueCount{T,1,UseGap})
	n.marginals[:] = n.counts
	n.total = sum(n.counts)
	n
end

function update!{T,UseGap}(n::ResidueCount{T,2,UseGap})
	n.marginals[:,1] = sum(n.counts, 2) # this is faster than sum(n,2)
	n.marginals[:,2] = sum(n.counts, 1)
	n.total = sum(n.marginals[:,1])
	n
end

function update!{T, N, UseGap}(n::ResidueCount{T, N, UseGap})
	for i in 1:N
		n.marginals[:,i] = sum(n.counts, _tuple_without_index(i,N))
	end
	n.total = sum(n.marginals[:,1])
	n
end

### Apply Pseudocount

# This is faster than array[:] += value
function __sum!(array, value)
  @inbounds for i in eachindex(array)
    array[i] += value
  end
  array
end

for (dim, gap, margi_exp, total_exp) in [ (:1, :true,  :(pse.λ), :(pse.λ * 21)),
																					(:1, :false, :(pse.λ), :(pse.λ * 20)),
																					(:2, :true,  :(pse.λ * 21), :(pse.λ * 441)),
																					(:2, :false, :(pse.λ * 20), :(pse.λ * 400)) ]
	@eval begin

		function apply_pseudocount!{T}(n::ResidueCount{T, $(dim), $(gap)}, pse::AdditiveSmoothing{T})
			margi_sum = $(margi_exp)
			total_sum = $(total_exp)
			__sum!(n.counts, pse.λ)
			__sum!(n.marginals, margi_sum)
			n.total += total_sum
			n
		end

		function fill!{T}(n::ResidueCount{T, $(dim), $(gap)}, pse::AdditiveSmoothing{T})
			margi_sum = $(margi_exp)
			total_sum = $(total_exp)
			fill!(n.counts, pse.λ)
			fill!(n.marginals, margi_sum)
			n.total = total_sum
			n
		end

	end

end

function apply_pseudocount!{T, N, UseGap}(n::ResidueCount{T, N, UseGap}, pse::AdditiveSmoothing{T})
	nres = UseGap ? 21 : 20
	margi_sum = pse.λ * (nres^(N-1))
	total_sum = pse.λ * (nres^N)
	__sum!(n.counts, pse.λ)
	__sum!(n.marginals, margi_sum)
	n.total += total_sum
	n
end

function fill!{T, N, UseGap}(n::ResidueCount{T, N, UseGap}, pse::AdditiveSmoothing{T})
	nres = UseGap ? 21 : 20
	margi_sum = pse.λ * (nres^(N-1))
	total_sum = pse.λ * (nres^N)
	fill!(n.counts, pse.λ)
	fill!(n.marginals, margi_sum)
	n.total = total_sum
	n
end

### Counting

for (usegap, testgap) in [ (:(true),:()), (:(false),:(aa == Int(GAP) && continue)) ]

  @eval begin

    function count!{T}(n::ResidueCount{T, 1, $(usegap)}, cl, res::AbstractVector{Residue})
      for i in 1:length(res)
        aa = Int(res[i])
        $(testgap)
        n.counts[aa] += getweight(cl,i)
      end
      update!(n)
    end

    count!{T}(n::ResidueCount{T, 1, $(usegap)}, res::AbstractVector{Residue}) = count!(n, one(T), res)

  end
end

function count!{T}(n::ResidueCount{T, 2, true}, cl, res1::AbstractVector{Residue}, res2::AbstractVector{Residue})
  for i in 1:length(res1)
    n.counts[Int(res1[i]), Int(res2[i])] += getweight(cl,i)
  end
  update!(n)
end

function count!{T}(n::ResidueCount{T, 2, false}, cl, res1::AbstractVector{Residue}, res2::AbstractVector{Residue})
  for i in 1:length(res1)
    aa1 = res1[i]
    aa2 = res2[i]
    if (aa1 != GAP) && (aa2 != GAP)
      n.counts[Int(aa1), Int(aa2)] += getweight(cl,i)
    end
  end
  update!(n)
end

count!{T}(n::ResidueCount{T, 2, true}, res1::AbstractVector{Residue}, res2::AbstractVector{Residue}) = count!(n, one(T), res1, res2)
count!{T}(n::ResidueCount{T, 2, false},res1::AbstractVector{Residue}, res2::AbstractVector{Residue}) = count!(n, one(T), res1, res2)

function count!{T, N, UseGap}(n::ResidueCount{T, N, UseGap}, cl, res::AbstractVector{Residue}...)
  if length(res) == N
    for i in 1:length(res[1])
      aa_list = [ Int(aa[i]) for aa in res ]
      if UseGap || (findfirst(aa_list, Int(GAP)) == 0)
        n.counts[aa_list...] += getweight(cl,i)
      end
    end
    update!(n)
  else
    throw("Number of arrays ($(length(res))) != $N")
  end
end

count!{T, N, UseGap}(n::ResidueCount{T, N, UseGap}, res::AbstractVector{Residue}...) = count!(n, one(T), res...)

#### Default counters

count(usegap::Bool, cl::Clusters, res::AbstractVector{Residue}...) = count!(zeros(ResidueCount{Float64, length(res), usegap}), cl, res...)
count(usegap::Bool, res::AbstractVector{Residue}...) = count!(zeros(ResidueCount{Int, length(res), usegap}), 1, res...)

## Probabilities

type ResidueProbability{N, UseGap} <: AbstractArray{Float64, N}
  probabilities::Array{Float64, N}
  marginals::Array{Float64, 2}
end

call(::Type{ResidueProbability{1, true}}) = ResidueProbability{1, true}(Array(Float64, 21), Array(Float64, (21,1)))
call(::Type{ResidueProbability{2, true}}) = ResidueProbability{2, true}(Array(Float64, (21, 21)), Array(Float64, (21,2)))

call(::Type{ResidueProbability{1, false}}) = ResidueProbability{1, false}(Array(Float64, 20), Array(Float64, (20,1)))
call(::Type{ResidueProbability{2, false}}) = ResidueProbability{2, false}(Array(Float64, (20, 20)), Array(Float64, (20,2)))

function call{N, UseGap}(::Type{ResidueProbability{N, UseGap}})
  nres = UseGap ? 21 : 20
  ResidueProbability{N, UseGap}(Array(Float64, (Int[ nres for d in 1:N]...)), Array(Float64, (nres, N)))
end

zeros(::Type{ResidueProbability{1, true}}) = ResidueProbability{1, true}(zeros(Float64, 21), zeros(Float64, (21, 1)))
zeros(::Type{ResidueProbability{2, true}}) = ResidueProbability{2, true}(zeros(Float64, (21, 21)), zeros(Float64, (21, 2)))

zeros(::Type{ResidueProbability{1, false}}) = ResidueProbability{1, false}(zeros(Float64, 20), zeros(Float64, (20, 1)))
zeros(::Type{ResidueProbability{2, false}}) = ResidueProbability{2, false}(zeros(Float64, (20, 20)), zeros(Float64, (20, 2)))

function zeros{N, UseGap}(::Type{ResidueProbability{N, UseGap}})
  nres = UseGap ? 21 : 20
  ResidueProbability{N, UseGap}(zeros(Float64, (Int[ nres for d in 1:N]...)), zeros(Float64, (nres, N)))
end

### Abstract Array Interface

Base.linearindexing(::ResidueProbability)	= Base.LinearFast()

"""`getindex(p::ResidueProbability, i::Int)` gives you access to the probabilities, use `Int()` for indexing with `Residues`"""
getindex(p::ResidueProbability, i::Int)	= getindex(p.probabilities, i)

"""`setindex!(p::ResidueProbability, v, i::Int)` set a value into probabilities field, but doesn't update the marginals.
Use `Int()` for indexing with `Residues`. Use `update!` in order to calculate again the marginals."""
setindex!(p::ResidueProbability, v, i::Int) = setindex!(p.probabilities, v, i)

#### Length & Size

length(p::ResidueProbability{1, true}) = 21
length(p::ResidueProbability{1,false}) = 20

length(p::ResidueProbability{2, true}) = 441
length(p::ResidueProbability{2,false}) = 400

length{N,UseGap}(p::ResidueProbability{N, UseGap}) = length(p.probabilities)

size(p::ResidueProbability{1, true}) = (21,)
size(p::ResidueProbability{1,false}) = (20,)

size(p::ResidueProbability{2, true}) = (21, 21)
size(p::ResidueProbability{2,false}) = (20, 20)

size{N,UseGap}(p::ResidueProbability{N, UseGap}) = size(p.probabilities)

#### Iteration Interface

start(p::ResidueProbability) = 1

@inbounds next(p::ResidueProbability, state::Int) = (p.probabilities[state], state + 1)

done{N,UseGap}(p::ResidueProbability{N, UseGap}, state) = state > length(p)

eltype{N,UseGap}(::Type{ResidueProbability{N, UseGap}}) = Float64

#### Similar

similar{N,UseGap}(p::ResidueProbability{N,UseGap}) = ResidueProbability{N,UseGap}()
similar{N,UseGap}(p::ResidueProbability{N,UseGap}, D::Int) = ResidueProbability{D,UseGap}()

### Update!

function update!{UseGap}(p::ResidueProbability{1,UseGap})
	p.marginals[:] = p.probabilities
	p
end

function update!{UseGap}(p::ResidueProbability{2,UseGap})
	p.marginals[:,1] = sum(p.probabilities, 2) # this is faster than sum(p,2)
	p.marginals[:,2] = sum(p.probabilities, 1)
	p
end

function update!{N, UseGap}(p::ResidueProbability{N, UseGap})
	for i in 1:N
		p.marginals[:,i] = sum(p.probabilities, _tuple_without_index(i,N))
	end
	p
end

### Normalize

# This is faster than array[:] /= value
function __div!(array, value)
  @inbounds for i in eachindex(array)
    array[i] /= value
  end
  array
end

"""```normalize!(p::ResidueProbability; updated::Bool=false)```

This function makes the sum of the probabilities to be one.
The sum is calculated using the `probabilities` field by default (It is assumed that the marginal are not updated).
The marginals are updated in the normalization.

If the marginals are updated, you can use `updated=true` for a faster normalization.
"""
function normalize!(p::ResidueProbability; updated::Bool=false)
  if !updated
    update!(p)
  end
  total = sum(p.marginals[:,1])
  if total != 1.0
    __div!(p.probabilities, total)
    __div!(p.marginals, total)
  end
  p
end

### ResidueProbability calculation from ResidueCount and others

# This is faster than p[:] = ( n ./ total )
function __fill_probabilities!{T,N}(p::Array{Float64,N}, n::Array{T,N}, total::T)
  @inbounds for i in 1:length(p) # p and n should have the same length
    p[i] = ( n[i] / total )
  end
  p
end

"""```fill!{T, N, UseGap}(p::ResidueProbability{N, UseGap}, n::ResidueCount{T, N, UseGap}; updated::Bool=false)```

This function fills a preallocated `ResidueProbability` (`p`) with the probabilities calculated from `n` (`ResidueCount`). This function updates `n` unless `updated=true`.

If `n` is updated, you can use `updated=true` for a faster calculation.
"""
function fill!{T, N, UseGap}(p::ResidueProbability{N, UseGap}, n::ResidueCount{T, N, UseGap}; updated::Bool=false)
  if !updated
    update!(n)
  end
  __fill_probabilities!(p.probabilities, n.counts, n.total)
  __fill_probabilities!(p.marginals, n.marginals, n.total)
	p
end

### BLOSUM based pseudofrequencies

function blosum_pseudofrequencies!(Gab::ResidueProbability{2,false}, Pab::ResidueProbability{2,false})
  @inbounds for a in 1:20, b in 1:20
    Gab[a,b] = zero(Float64)
    for i in 1:20, j in 1:20
      P = Pab[i,j]
      if P != 0
	      # BLOSUM62_P_i_j[i,a] is p(a | i)
	      Gab[a,b] += ( P * BLOSUM62_Pij[i,a] * BLOSUM62_Pij[j,b] )
      end
    end
  end
  normalize!(Gab)
end

### Apply pseudofrequencies

function apply_pseudofrequencies!(Pab::ResidueProbability{2,false}, Gab::ResidueProbability{2,false}, α, β)
  frac = one(Float64) / ( α + β )
  total = zero(Float64)
  @inbounds for i in 1:20, j in 1:20
    value = ( α * Pab[i,j] + β * Gab[i,j] ) * frac
    Pab[i,j] = value
    total += value
  end
  if total != 1.0
    __div!(Pab.probabilities, total)
    update!(Pab)
  else
    update!(Pab)
  end
end

# """```probabilities{T, N, UseGap}(n::ResidueCount{T, N, UseGap})```

# This function creates a new `ResidueProbability` matrix based on the counts of `n`. If you want to calculate probabilities on a preallocated `ResidueProbability` matrix use `fill!` instead.
# """
# function probabilities{T, N, UseGap}(n::ResidueCount{T, N, UseGap})
#   p = ResidueProbability{N, UseGap}()
#   __fill_probabilities!(p.probabilities, n.counts, sum(n.counts)) # slow but safe: n.total will be wrong if ResidueCount is not updated
#   update!(p)
# end

# # Use it only when n::ResidueCount is updated
# function __probabilities{T, N, UseGap}(n::ResidueCount{T, N, UseGap})
#   p = ResidueProbability{N, UseGap}()
#   __fill_probabilities!(p.probabilities, n.counts, n.total)
#   __fill_probabilities!(p.marginals, n.marginals, n.total)
# 	p
# end
# """```probabilities{T, N, UseGap}(usegap::Bool[, cl::Clusters], res::AbstractVector{Residue}...)```

# This function creates a new `ResidueProbability` matrix based on the counts of residues in the sequences/positions `res`.
# i.e. If you use two columns of a MSA, you obtain the probabilities for pairs of residues between the two columns """
# probabilities(usegap::Bool, cl::Clusters, res::AbstractVector{Residue}...) = __probabilities(count(usegap, cl, res...))
# probabilities(usegap::Bool, res::AbstractVector{Residue}...) = __probabilities(count(usegap, res...))
