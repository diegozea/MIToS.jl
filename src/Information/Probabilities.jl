import Base: zero, zeros, start, next, done, length, eltype, 
						 size, setindex!, 
						 getindex, similar #, print # , copy, deepcopy, fill!, getindex, setindex!

# Counts and Pseudocounts

abstract Pseudocount

immutable AdditiveSmoothing{T} <: Pseudocount
  λ::T
end

zero{T}(::Type{AdditiveSmoothing{T}}) = AdditiveSmoothing{T}(zero(T))

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

call{T}(::Type{ResidueCount{T, 1, true}}) = ResidueCount{T, 1, true}(Array(T, 21), Array(T, (1,21)), zero(T))
call{T}(::Type{ResidueCount{T, 2, true}}) = ResidueCount{T, 2, true}(Array(T, (21, 21)), Array(T, (2,21)), zero(T))

call{T}(::Type{ResidueCount{T, 1, false}}) = ResidueCount{T, 1, false}(Array(T, 20), Array(T, (1,20)), zero(T))
call{T}(::Type{ResidueCount{T, 2, false}}) = ResidueCount{T, 2, false}(Array(T, (20, 20)), Array(T, (2,20)), zero(T))

function call{T, N, UseGap}(::Type{ResidueCount{T, N, UseGap}})
  n = UseGap ? 21 : 20
  ResidueCount{T, N, UseGap}(Array(T, (Int[ n for d in 1:N]...)), Array(T, (N, n)), zero(T))
end

zeros{T}(::Type{ResidueCount{T, 1, true}}) = ResidueCount{T, 1, true}(zeros(T, 21), zeros(T, (1,21)), zero(T))
zeros{T}(::Type{ResidueCount{T, 2, true}}) = ResidueCount{T, 2, true}(zeros(T, (21, 21)), zeros(T, (2,21)), zero(T))

zeros{T}(::Type{ResidueCount{T, 1, false}}) = ResidueCount{T, 1, false}(zeros(T, 20), zeros(T, (1,20)), zero(T))
zeros{T}(::Type{ResidueCount{T, 2, false}}) = ResidueCount{T, 2, false}(zeros(T, (20, 20)), zeros(T, (2,20)), zero(T))

function zeros{T, N, UseGap}(::Type{ResidueCount{T, N, UseGap}})
  n = UseGap ? 21 : 20
  ResidueCount{T, N, UseGap}(zeros(T, (Int[ n for d in 1:N]...)), zeros(T, (N, n)), zero(T))
end

### Abstract Array Interface

Base.linearindexing(::ResidueCount)	= Base.LinearFast()

"""`getindex(n::ResidueCount, i::Int)` gives you access to the counts, use `Int()` for indexing with `Residues`"""
getindex(n::ResidueCount, i::Int)	= getindex(n.counts, i)

"""`setindex!(n::ResidueCount, v, i::Int)` set a value into counts field, but doesn't update the marginals and total.
Use `Int()` for indexing with `Residues`. Use `update!` in order to calculate again the marginals and total."""
setindex!(n::ResidueCount, v, i::Int) = setindex!(N.counts, v, i)

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


# ### PROBABILITIES AND PSEUDOCOUNTS ###
# 
# abstract Pseudocount
# 
# immutable Fixed <: Pseudocount
#   λ::Float64
# end
# 
# type Pseudofrequencies <: Pseudocount
#   α::Float64
#   β::Float64
#   Gab::Matrix{Float64}
# end
# 
# type ResidueProbabilities
#   Pa::Vector{Float64}
# end
# 
# type ResiduePairProbabilities
#   Pab::Matrix{Float64}
#   Pa::Vector{Float64}
#   Pb::Vector{Float64}
# end
# 
# ### Constructors ###
# 
# zero(::Type{Fixed}) = Fixed(zero(Float64))
# 
# zeros(::Type{Pseudofrequencies}) = Pseudofrequencies(zero(Float64), zero(Float64), zeros(Float64, (20,20)))
# zeros(::Type{ResidueProbabilities}) = ResidueProbabilities(zeros(Float64,20))
# zeros(::Type{ResiduePairProbabilities}) = ResiduePairProbabilities(zeros(Float64, (20,20)), zeros(Float64,20), zeros(Float64,20))
# 
# Fixed() = zero(Fixed)
# Pseudofrequencies() = Pseudofrequencies(zero(Float64), zero(Float64), Array(Float64, (20,20)))
# ResidueProbabilities() = ResidueProbabilities(Array(Float64,20))
# ResiduePairProbabilities() = ResiduePairProbabilities(Array(Float64, (20,20)), Array(Float64,20), Array(Float64,20))
# 
# Pseudofrequencies(α::Float64, β::Float64) = Pseudofrequencies(α, β, Array(Float64, (20,20)))
# 
# ### Copy ###
# 
# copy(x::Fixed) = Fixed(copy(x.λ))
# copy(x::Pseudofrequencies) = Pseudofrequencies(copy(x.α), copy(x.β), copy(x.Gab))
# copy(pa::ResidueProbabilities) = ResidueProbabilities(copy(pa.Pa))
# copy(pab::ResiduePairProbabilities) = ResiduePairProbabilities(copy(pab.Pab), copy(pab.Pa), copy(pab.Pb))
# 
# ### Indexing ###
# 
# getindex(pa::ResidueProbabilities, i::Int) = getindex(pa.Pa, i)
# setindex!(pa::ResidueProbabilities, p::Float64, i::Int) = setindex!(pa.Pa, p, i)
# 
# getindex(pab::ResiduePairProbabilities, i::Int, j::Int) = getindex(pab.Pab, i, j)
# 
# getindex(pse::Pseudofrequencies, i::Int, j::Int) = getindex(pse.Gab, i, j)
# 
# ### Estimation (Filling the Probabilities) ###
# 
# """Fill the Pab matrix with λ, but doesn't update the marginal probabilities"""
# function __initialize!(pab::ResiduePairProbabilities)
#   fill!(pab.Pab, zero(Float64))
#   zero(Float64)
# end
# 
# function __initialize!(pab::ResiduePairProbabilities, pseudocount::Fixed)
#   fill!(pab.Pab, pseudocount.λ)
#   pseudocount.λ * 400.0
# end
# 
# function __finalize!(pab::ResiduePairProbabilities, total::Float64)
#   pab.Pab[:,:] /= total
#   pab.Pa[:] = sum(pab.Pab,1)
#   pab.Pa[:] = sum(pab.Pab,2)
#   pab
# end
# 
# function __fill_kernel!(pab::ResiduePairProbabilities, seqi::AbstractVector{Residue}, seqj::AbstractVector{Residue}, total::Float64)
#   for i in 1:length(seqi) # seqj should have the same length
#     if seqi[i] != GAP  && seqj[i] != GAP
#       @inbounds pab.Pab[ seqi[i] , seqj[i] ] += one(Float64)
#       total += one(Float64)
#     end
#   end
#   total
# end
# 
# function __fill_kernel!(pab::ResiduePairProbabilities, cl::Clusters, seqi::AbstractVector{Residue}, seqj::AbstractVector{Residue}, total::Float64)
#   for i in 1:length(seqi) # seqj should have the same length
#     if seqi[i] != GAP  && seqj[i] != GAP
#       weight = getweight(cl, i)
#       @inbounds pab.Pab[ seqi[i] , seqj[i] ] += weight
#       total += weight
#     end
#   end
#   total
# end
# 
# function __fill_pseudofrequencies!(pseudocount::Pseudofrequencies, pab::ResiduePairProbabilities)
#   total = zero(Float64)
#   for a in 1:20, b in 1:20
#     @inbounds pseudocount.Gab[a,b] = zero(Float64)
#     for i in 1:20, j in 1:20
#       P = pab[i,j]
#       if P != 0
# 	      # BLOSUM62_P_i_j[i,a] is p(a | i)
# 	      @inbounds val = P * BLOSUM62_Pij[i,a] * BLOSUM62_Pij[j,b]
# 	      @inbounds pseudocount.Gab[a,b] += val
# 	      total += val
#       end
#     end
#   end
#   pseudocount.Gab[:,:] /= total
#   pseudocount
# end
# 
# function __apply_pseudofrequencies!(pab::ResiduePairProbabilities, pseudocount::Pseudofrequencies)
#   frac = one(Float64) / ( pseudocount.α + pseudocount.β )
#   total = zero(Float64)
#   for i in 1:20, j in 1:20
#     @inbounds val = ( pseudocount.α * pab[i,j] + pseudocount.β * pseudocount[i,j] ) * frac
#     @inbounds pab.Pab[i,j] = val
#     total += val
#   end
#   total
# end
# 
# function fill!(pab::ResiduePairProbabilities, seqi::AbstractVector{Residue}, seqj::AbstractVector{Residue})
#   total = __initialize!(pab)
#   total = __fill_kernel!(pab, seqi, seqj, total)
#   __finalize!(pab, total)
# end
# 
# function fill!(pab::ResiduePairProbabilities, pseudocount::Fixed, seqi::AbstractVector{Residue}, seqj::AbstractVector{Residue})
#   total = __initialize!(pab, pseudocount)
#   total = __fill_kernel!(pab, seqi, seqj, total)
#   __finalize!(pab, total)
# end
# 
# function fill!(pab::ResiduePairProbabilities, pseudocount::Fixed, cl::Clusters, seqi::AbstractVector{Residue}, seqj::AbstractVector{Residue})
#   total = __initialize!(pab, pseudocount)
#   total = __fill_kernel!(pab, cl, seqi, seqj, total)
#   __finalize!(pab, total)
# end
# 
# function fill!(pab::ResiduePairProbabilities, pseudocount::Pseudofrequencies, seqi::AbstractVector{Residue}, seqj::AbstractVector{Residue})
#   total = __initialize!(pab)
#   total = __fill_kernel!(pab, seqi, seqj, total)
#   if pseudocount.β != zero(Float64)
#     __fill_pseudofrequencies!(pseudocount, pab)
#     total = __apply_pseudofrequencies(pab, pseudocount)
#   end
#   __finalize!(pab, total)
# end
# 
# function fill!(pab::ResiduePairProbabilities, pseudocount::Pseudofrequencies, cl::Clusters, seqi::AbstractVector{Residue}, seqj::AbstractVector{Residue})
#   total = __initialize!(pab)
#   total = __fill_kernel!(pab, cl, seqi, seqj, total)
#   if pseudocount.β != zero(Float64)
#     __fill_pseudofrequencies!(pseudocount, pab)
#     total = __apply_pseudofrequencies(pab, pseudocount)
#   end
#   __finalize!(pab, total)
# end
