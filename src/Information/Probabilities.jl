"""
`SequenceWeights` is an alias for `Union{ClusteringResult, AbstractVector}`.

This is the type of the keyword argment `weight` of `count`, `probabilities` and derived functions.
The type should define `getweight` to be useful in those functions.
"""
typealias SequenceWeights Union{ClusteringResult, AbstractVector} # AbstractVector because getweight(cl) returns a Vector

"""
Number of residues used in the `ResidueContingencyTables`.
A 20x20 table returns 20 (ungapped alphabet).
"""
@inline nresidues{T, N}(n::ResidueContingencyTables{T, N, true})  = 21
@inline nresidues{T, N}(n::ResidueContingencyTables{T, N, false}) = 20

## Counts

"""
`ResidueCount{T, N, UseGap}` is used to count residues in columns (or sequences) of an MSA.
`N` is the dimensionality and should be an `Int`, i.e. 2 if 2 columns are used to count pairs.
`UseGap` is a `Bool`, `true` means that **ResidueCount** counts gaps in the position 21.

- The field marginal is used for pre allocation of marginal sums.
- The field total is used for storing the table sum.
"""
type ResidueCount{T, N, UseGap} <: ResidueContingencyTables{T, N, UseGap}
  table::Array{T, N}
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

#### Similar

similar{T,N,UseGap}(n::ResidueCount{T,N,UseGap}) = ResidueCount{T,N,UseGap}()
similar{T,N,S,UseGap}(n::ResidueCount{T,N,UseGap}, ::Type{S}) = ResidueCount{S,N,UseGap}()
similar{T,N,UseGap}(n::ResidueCount{T,N,UseGap}, D::Int) = ResidueCount{T,D,UseGap}()
similar{T,N,S,UseGap}(n::ResidueCount{T,N,UseGap}, ::Type{S}, D::Int) = ResidueCount{S,D,UseGap}()

### Update!

"Updates the marginals values (and total value for `ResidueCount`) of a `ResidueContingencyTables`"
function update!{T,UseGap}(n::ResidueCount{T,1,UseGap})
	n.marginals[:] = n.table
	n.total = sum(n.table)
	n
end

function update!{T,UseGap}(n::ResidueCount{T,2,UseGap})
	n.marginals[:,1] = sum(n.table, 2) # this is faster than sum(n,2)
	n.marginals[:,2] = sum(n.table, 1)
	n.total = sum(n.marginals[:,1])
	n
end

function update!{T, N, UseGap}(n::ResidueCount{T, N, UseGap})
	@inbounds for i in 1:N
		n.marginals[:,i] = sum(n.table, _tuple_without_index(i,N))
	end
	n.total = sum(n.marginals[:,1])
	n
end

### Apply Pseudocount

# This is faster than array[:] += value

"""
`apply_pseudocount!{T, N, UseGap}(n::ResidueCount{T, N, UseGap}, pse::AdditiveSmoothing{T})`

Uses an instance of `AdditiveSmoothing` to efficiently fill with a constant value each element of the table.
"""
function apply_pseudocount!{T, N, UseGap}(n::ResidueCount{T, N, UseGap}, pse::AdditiveSmoothing{T})
	nres = nresidues(n)
	margi_sum = pse.λ * (nres^(N-1))
	total_sum = pse.λ * (nres^N)
	_sum!(n.table, pse.λ)
	_sum!(n.marginals, margi_sum)
	n.total += total_sum
	n
end

"""
`fill!{T, N, UseGap}(n::ResidueCount{T, N, UseGap}, pse::AdditiveSmoothing{T})`
fills a preallocated `ResidueCount` (`p`) with the pseudocount (`pse`).
"""
function fill!{T, N, UseGap}(n::ResidueCount{T, N, UseGap}, pse::AdditiveSmoothing{T})
	nres = nresidues(n)
	margi_sum = pse.λ * (nres^(N-1))
	total_sum = pse.λ * (nres^N)
	fill!(n.table, pse.λ)
	fill!(n.marginals, margi_sum)
	n.total = total_sum
	n
end

### Counting

"""
`count!` adds counts from vector of residues to a `ResidueCount` object.
It can take a SequenceWeights object as second argument.

```julia

julia> using MIToS.Information

julia> using MIToS.MSA

julia> seq = Residue[ i for i in 1:20];

julia> Ni = count(seq)
20-element MIToS.Information.ResidueCount{Float64,1,false}:
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0

julia> count!(Ni, seq)
20-element MIToS.Information.ResidueCount{Float64,1,false}:
 2.0
 2.0
 2.0
 2.0
 2.0
 2.0
 2.0
 2.0
 2.0
 2.0
 2.0
 2.0
 2.0
 2.0
 2.0
 2.0
 2.0
 2.0
 2.0
 2.0

```
"""
count!{T, N, UseGap}(n::ResidueCount{T, N, UseGap}, res::AbstractVector{Residue}...) = count!(n, NoClustering(), res...)

#### Default counters

"""
```
count(res::AbstractVector{Residue}...; usegap=false, weight=NoClustering())
count(pseudocount::Pseudocount, res::AbstractVector{Residue}...; usegap=false, weight=NoClustering())
```

`count` creates a new ResidueCount counting the number of residues, pairs of residues, etc. in the sequences/columns.
"""
count(res::AbstractVector{Residue}...; usegap::Bool=false,
	weight::SequenceWeights=NoClustering()) = count!(zeros(ResidueCount{Float64, length(res), usegap}), weight, res...)

count(pseudocount::Pseudocount, res::AbstractVector{Residue}...; usegap::Bool=false,
	weight::SequenceWeights=NoClustering()) = apply_pseudocount!(count!(zeros(ResidueCount{Float64, length(res), usegap}), weight, res...), pseudocount)

count{T}(pseudocount::AdditiveSmoothing{T}, res::AbstractVector{Residue}...; usegap::Bool=false,
	weight::SequenceWeights=NoClustering()) = apply_pseudocount!(count!(zeros(ResidueCount{T, length(res), usegap}), weight, res...), pseudocount)

## Probabilities

"""
`ResidueProbability{T, N, UseGap}` is used to store residue probabilities.
`N` is the dimensionality and should be an `Int`, e.g. 2 to store probabilities of residue pairs.
`UseGap` is a `Bool`, `true` means that gap probabilities are stored in the position 21 of each dimension.

- The field marginal is used for pre allocation of marginal sums.
"""
type ResidueProbability{T, N, UseGap} <: ResidueContingencyTables{T, N, UseGap}
  table::Array{T, N}
  marginals::Array{T, 2}
end

call{T}(::Type{ResidueProbability{T, 1, true}}) = ResidueProbability{T, 1, true}(Array(T, 21), Array(T, (21,1)))
call{T}(::Type{ResidueProbability{T, 2, true}}) = ResidueProbability{T, 2, true}(Array(T, (21, 21)), Array(T, (21,2)))

call{T}(::Type{ResidueProbability{T, 1, false}}) = ResidueProbability{T, 1, false}(Array(T, 20), Array(T, (20,1)))
call{T}(::Type{ResidueProbability{T, 2, false}}) = ResidueProbability{T, 2, false}(Array(T, (20, 20)), Array(T, (20,2)))

function call{T, N, UseGap}(::Type{ResidueProbability{T, N, UseGap}})
  nres = UseGap ? 21 : 20
  ResidueProbability{T, N, UseGap}(Array(T, (Int[ nres for d in 1:N]...)), Array(T, (nres, N)))
end

zeros{T}(::Type{ResidueProbability{T, 1, true}}) = ResidueProbability{T, 1, true}(zeros(T, 21), zeros(T, (21, 1)))
zeros{T}(::Type{ResidueProbability{T, 2, true}}) = ResidueProbability{T, 2, true}(zeros(T, (21, 21)), zeros(T, (21, 2)))

zeros{T}(::Type{ResidueProbability{T, 1, false}}) = ResidueProbability{T, 1, false}(zeros(T, 20), zeros(T, (20, 1)))
zeros{T}(::Type{ResidueProbability{T, 2, false}}) = ResidueProbability{T, 2, false}(zeros(T, (20, 20)), zeros(T, (20, 2)))

function zeros{T, N, UseGap}(::Type{ResidueProbability{T, N, UseGap}})
  nres = UseGap ? 21 : 20
  ResidueProbability{T, N, UseGap}(zeros(T, (Int[ nres for d in 1:N]...)), zeros(T, (nres, N)))
end

#### Similar

similar{T, N,UseGap}(p::ResidueProbability{T, N,UseGap}) = ResidueProbability{T, N,UseGap}()
similar{T, N,UseGap}(p::ResidueProbability{T, N,UseGap}, D::Int) = ResidueProbability{T, D,UseGap}()

## Default probabilities

"""
```
probabilities(res::AbstractVector{Residue}...; usegap=false, weight=NoClustering())
probabilities(pseudocount::Pseudocount, res::AbstractVector{Residue}...; usegap=false, weight=NoClustering())
```

`probabilities` creates a new ResidueProbability with the probabilities of residues, pairs of residues, etc. in the sequences/columns.
"""
probabilities{T}(::Type{T}, res::AbstractVector{Residue}...; usegap::Bool=false,
	weight::SequenceWeights=NoClustering()) = fill!(ResidueProbability{T, length(res), usegap}(), count(res..., usegap=usegap, weight=weight))

probabilities{T}(::Type{T}, pseudocount::Pseudocount, res::AbstractVector{Residue}...; usegap::Bool=false,
	weight::SequenceWeights=NoClustering()) = fill!(ResidueProbability{T, length(res), usegap}(), count(pseudocount, res..., usegap=usegap, weight=weight))

"""
`probabilities(T, α, β, res1, res2, [weight])` use BLOSUM62 based pseudofrequencies.
α is the weight of the evidence, and β the weight of the pseudofrequencies.
"""
function probabilities{T}(::Type{T}, α, β, res1::AbstractVector{Residue}, res2::AbstractVector{Residue}; weight::SequenceWeights=NoClustering())
	Pab = fill!(ResidueProbability{T, 2, false}(), count(res1, res2, usegap=false, weight=weight))
	Gab = blosum_pseudofrequencies!(ResidueProbability{T, 2,false}(), Pab)
	apply_pseudofrequencies!(Pab, Gab, α, β)
end
