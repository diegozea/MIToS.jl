import  Base: zero, one, zeros, start, next, done, length, eltype,
        size, setindex!, getindex, similar, fill!, count #, print # , copy, deepcopy, fill!, getindex, setindex!

"""
`SequenceWeights` is an alias for `Union{ClusteringResult, AbstractVector}`.

This is the type of the keyword argment `weight` of `count`, `probabilities` and derived functions.
The type should define `getweight` to be useful in those functions.
"""
typealias SequenceWeights Union{ClusteringResult, AbstractVector} # AbstractVector because getweight(cl) returns a Vector


## ResidueContingencyTables

abstract ResidueContingencyTables{T, N, UseGap} <: AbstractArray{T, N}

### Abstract Array Interface

Base.linearindexing(::ResidueContingencyTables)	= Base.LinearFast()

"""
`getindex(n::ResidueContingencyTables, i::Int)`

Gives access to the table. Use `Int()` for indexing with `Residue`s.
"""
getindex(n::ResidueContingencyTables, i::Int)	= getindex(n.table, i)

"""
`setindex!(n::ResidueCount, v, i::Int)` set a value into the table, but doesn't update the marginals (and total).
Use `Int()` for indexing with `Residues`. Use `update!` in order to calculate again the marginals (and total).
"""
setindex!(n::ResidueContingencyTables, v, i::Int) = setindex!(n.table, v, i)

#### Length & Size

length{T}(n::ResidueContingencyTables{T, 1, true}) = 21
length{T}(n::ResidueContingencyTables{T, 1,false}) = 20

length{T}(n::ResidueContingencyTables{T, 2, true}) = 441
length{T}(n::ResidueContingencyTables{T, 2,false}) = 400

length{T,N,UseGap}(n::ResidueContingencyTables{T, N, UseGap}) = length(n.table)

size{T}(n::ResidueContingencyTables{T, 1, true}) = (21,)
size{T}(n::ResidueContingencyTables{T, 1,false}) = (20,)

size{T}(n::ResidueContingencyTables{T, 2, true}) = (21, 21)
size{T}(n::ResidueContingencyTables{T, 2,false}) = (20, 20)

size{T,N,UseGap}(n::ResidueContingencyTables{T, N, UseGap}) = size(n.table)

"""
Number of residues used in the `ResidueContingencyTables`.
A 20x20 table returns 20 (ungapped alphabet).
"""
@inline nresidues{T, N}(n::ResidueContingencyTables{T, N, true})  = 21
@inline nresidues{T, N}(n::ResidueContingencyTables{T, N, false}) = 20

#### Iteration Interface

start(n::ResidueContingencyTables) = 1

@inbounds next(n::ResidueContingencyTables, state::Int) = (n.table[state], state + 1)

done{T,N,UseGap}(n::ResidueContingencyTables{T, N, UseGap}, state) = state > length(n)

eltype{T,N,UseGap}(::Type{ResidueContingencyTables{T, N, UseGap}}) = T

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
function _sum!(array, value)
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
			_sum!(n.table, pse.λ)
			_sum!(n.marginals, margi_sum)
			n.total += total_sum
			n
		end

		function fill!{T}(n::ResidueCount{T, $(dim), $(gap)}, pse::AdditiveSmoothing{T})
			margi_sum = $(margi_exp)
			total_sum = $(total_exp)
			fill!(n.table, pse.λ)
			fill!(n.marginals, margi_sum)
			n.total = total_sum
			n
		end

	end

end

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

function count!{T}(n::ResidueCount{T, 1, true}, weights, res::AbstractVector{Residue})
    @simd for i in 1:length(res)
        @inbounds n.table[Int(res[i])] += getweight(weights,i)
    end
    update!(n)
end

function count!{T}(n::ResidueCount{T, 1, false}, weights, res::AbstractVector{Residue})
    @inbounds for i in 1:length(res)
        aa = Int(res[i])
        aa == Int(GAP) && continue
        n.table[aa] += getweight(weights,i)
    end
    update!(n)
end

count!{T, UseGap}(n::ResidueCount{T, 1, UseGap}, res::AbstractVector{Residue}) = count!(n, NoClustering(), res)

function count!{T}(n::ResidueCount{T, 2, true}, weights, res1::AbstractVector{Residue}, res2::AbstractVector{Residue})
  @assert length(res1) == length(res2)
  @simd for i in 1:length(res1)
    @inbounds n.table[Int(res1[i]), Int(res2[i])] += getweight(weights,i)
  end
  update!(n)
end

function count!{T}(n::ResidueCount{T, 2, false}, weights, res1::AbstractVector{Residue}, res2::AbstractVector{Residue})
  for i in 1:length(res1)
    aa1 = res1[i]
    aa2 = res2[i]
    if (aa1 != GAP) && (aa2 != GAP)
      n.table[Int(aa1), Int(aa2)] += getweight(weights,i)
    end
  end
  update!(n)
end

count!{T}(n::ResidueCount{T, 2, true}, res1::AbstractVector{Residue}, res2::AbstractVector{Residue}) = count!(n, NoClustering(), res1, res2)
count!{T}(n::ResidueCount{T, 2, false},res1::AbstractVector{Residue}, res2::AbstractVector{Residue}) = count!(n, NoClustering(), res1, res2)

function count!{T, N, UseGap}(n::ResidueCount{T, N, UseGap}, weights, res::AbstractVector{Residue}...)
  if length(res) == N
    for i in 1:length(res[1])
      aa_list = [ Int(aa[i]) for aa in res ]
      if UseGap || (findfirst(aa_list, Int(GAP)) == 0)
        n.table[aa_list...] += getweight(weights,i)
      end
    end
    update!(n)
  else
    throw("Number of arrays ($(length(res))) != $N")
  end
end

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
