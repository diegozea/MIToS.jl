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
