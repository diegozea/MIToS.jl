function _get_matrix_residue(msa::AbstractMultipleSequenceAlignment)::Matrix{Residue}
    getarray(namedmatrix(msa))
end
_get_matrix_residue(msa::NamedArray)::Matrix{Residue} = getarray(msa)
_get_matrix_residue(msa::Matrix{Residue})::Matrix{Residue} = msa

# Kernel function to fill ContingencyTable based on residues
function _mapfreq_kernel!(
    f::Function,
    freqtable::Union{Probabilities{T,N,A},Frequencies{T,N,A}},
    weights,
    pseudocounts,
    pseudofrequencies,
    res::Vararg{AbstractVector{Residue},N};
    kargs...,
) where {T,N,A}
    table = freqtable.table
    _cleanup_table!(table) # frequencies! calls _cleanup_temporal! and cleans marginals
    frequencies!(table, res..., weights = weights, pseudocounts = pseudocounts) # frequencies! calls apply_pseudocount! and  _update!
    if isa(freqtable, Probabilities{T,N,A})
        normalize!(table)
        apply_pseudofrequencies!(table, pseudofrequencies)
    end
    f(freqtable; kargs...)
end

_mapfreq_kargs_doc = """
- `weights` (default: `NoClustering()`): Weights to be used for table counting.
- `pseudocounts` (default: `NoPseudocount()`): `Pseudocount` object to be applied to table.
- `pseudofrequencies` (default: `NoPseudofrequencies()`): `Pseudofrequencies` to be applied to the normalized (probabilities) table.
"""

_mappairfreq_kargs_doc = """
- `usediagonal` (default: `true`): If `true`, the function will be also applied to the diagonal elements.
- `diagonalvalue` (default: `zero`): Value to fill diagonal elements if `usediagonal` is `false`.
"""

# Residues: The output is a Named Vector
# --------------------------------------

function _mapfreq!(
    f::Function,
    res_list::Vector{V}, # sequences or columns
    table::Union{Probabilities{T,1,A},Frequencies{T,1,A}};
    weights::WeightTypes = NoClustering(),
    pseudocounts::Pseudocount = NoPseudocount(),
    pseudofrequencies::Pseudofrequencies = NoPseudofrequencies(),
    kargs...,
) where {T,A,V<:AbstractArray{Residue}}
    scores = map(res_list) do res
        _mapfreq_kernel!(f, table, weights, pseudocounts, pseudofrequencies, res; kargs...)
    end
    scores
end

# Map to each column

"""
It efficiently map a function (first argument) that takes a table of `Frequencies` or
`Probabilities` (third argument). The table is filled in place with the counts or
probabilities of each column from the `msa` (second argument).

$_mapfreq_kargs_doc
"""
function mapcolfreq!(
    f::Function,
    msa::AbstractMatrix{Residue},
    table::Union{Probabilities{T,1,A},Frequencies{T,1,A}};
    weights::WeightTypes = NoClustering(),
    pseudocounts::Pseudocount = NoPseudocount(),
    pseudofrequencies::Pseudofrequencies = NoPseudofrequencies(),
    kargs...,
) where {T,A}
    N = ncolumns(msa)
    residues = _get_matrix_residue(msa)
    ncol = size(residues, 2)
    columns = map(i -> view(residues, :, i), 1:ncol) # 2x faster than calling view inside the loop
    scores = _mapfreq!(
        f,
        columns,
        table;
        weights = weights,
        pseudocounts = pseudocounts,
        pseudofrequencies = pseudofrequencies,
        kargs...,
    )
    name_list = columnnames(msa)
    NamedArray(
        reshape(scores, (1, N)),
        (
            OrderedDict{String,Int}(string(f) => 1),
            OrderedDict{String,Int}(name_list[i] => i for i = 1:N),
        ),
        ("Function", "Col"),
    )
end

# Map to each sequence

"""
It efficiently map a function (first argument) that takes a table of `Frequencies` or
`Probabilities` (third argument). The table is filled in place with the counts or
probabilities of each sequence from the `msa` (second argument).

$_mapfreq_kargs_doc
"""
function mapseqfreq!(
    f::Function,
    msa::AbstractMatrix{Residue},
    table::Union{Probabilities{T,1,A},Frequencies{T,1,A}};
    weights::WeightTypes = NoClustering(),
    pseudocounts::Pseudocount = NoPseudocount(),
    pseudofrequencies::Pseudofrequencies = NoPseudofrequencies(),
    kargs...,
) where {T,A}
    N = nsequences(msa)
    sequences = getresiduesequences(msa)
    name_list = sequencenames(msa)
    scores = _mapfreq!(
        f,
        sequences,
        table;
        weights = weights,
        pseudocounts = pseudocounts,
        pseudofrequencies = pseudofrequencies,
        kargs...,
    )
    NamedArray(
        reshape(scores, (N, 1)),
        (
            OrderedDict{String,Int}(name_list[i] => i for i = 1:N),
            OrderedDict{String,Int}(string(f) => 1),
        ),
        ("Seq", "Function"),
    )
end

# Residue pairs: The output is a Named PairwiseListMatrix
# -------------------------------------------------------

function _mappairfreq!(
    f::Function,
    res_list::Vector{V}, # sequences or columns
    plm::PairwiseListMatrix{T,D,TV}, # output
    table::Union{Probabilities{T,2,A},Frequencies{T,2,A}},
    usediagonal::Type{Val{D}} = Val{true};
    weights::WeightTypes = NoClustering(),
    pseudocounts::Pseudocount = NoPseudocount(),
    pseudofrequencies::Pseudofrequencies = NoPseudofrequencies(),
    kargs...,
) where {T,D,TV,A,V<:AbstractArray{Residue}}
    @inbounds @iterateupper plm D begin
        list[k] = _mapfreq_kernel!(
            f,
            table,
            weights,
            pseudocounts,
            pseudofrequencies,
            res_list[i],
            res_list[j];
            kargs...,
        )
    end
    plm
end

# Map to column pairs

"""
It efficiently map a function (first argument) that takes a table of `Frequencies` or
`Probabilities` (third argument). The table is filled in place with the counts or
probabilities of each pair of columns from the `msa` (second argument).

$_mapfreq_kargs_doc
$_mappairfreq_kargs_doc
"""
function mapcolpairfreq!(
    f::Function,
    msa::AbstractMatrix{Residue},
    table::Union{Probabilities{T,2,A},Frequencies{T,2,A}};
    usediagonal::Bool = true,
    diagonalvalue::T = zero(T),
    weights::WeightTypes = NoClustering(),
    pseudocounts::Pseudocount = NoPseudocount(),
    pseudofrequencies::Pseudofrequencies = NoPseudofrequencies(),
    kargs...,
) where {T,A}
    ncol = ncolumns(msa)
    residues = _get_matrix_residue(msa)
    columns = map(i -> view(residues, :, i), 1:ncol) # 2x faster than calling view inside the loop
    scores = columnpairsmatrix(msa, T, Val{usediagonal}, diagonalvalue) # Named PairwiseListMatrix
    plm = getarray(scores)
    _mappairfreq!(
        f,
        columns,
        plm,
        table,
        Val{usediagonal};
        weights = weights,
        pseudocounts = pseudocounts,
        pseudofrequencies = pseudofrequencies,
        kargs...,
    )
    scores
end

# DEPRECATED, usediagonal is now a boolean keyword argument
function mapcolpairfreq!(f, msa, table, usediagonal::Type{Val{D}}; kargs...) where {D}
    Base.depwarn(
        "The `usediagonal` positional argument of `mapcolpairfreq!` taking `Val{true}` or `Val{false}` is deprecated. Use `usediagonal = true` or `usediagonal = false` instead.",
        :mapcolpairfreq!,
        force = true,
    )
    mapcolpairfreq!(f, msa, table; usediagonal = D, kargs...)
end

# Map to sequence pairs

"""
It efficiently map a function (first argument) that takes a table of `Frequencies` or
`Probabilities` (third argument). The table is filled in place with the counts or
probabilities of each pair of sequences from the `msa` (second argument).

$_mapfreq_kargs_doc
$_mappairfreq_kargs_doc
"""
function mapseqpairfreq!(
    f::Function,
    msa::AbstractMatrix{Residue},
    table::Union{Probabilities{T,2,A},Frequencies{T,2,A}};
    usediagonal::Bool = true,
    diagonalvalue::T = zero(T),
    weights::WeightTypes = NoClustering(),
    pseudocounts::Pseudocount = NoPseudocount(),
    pseudofrequencies::Pseudofrequencies = NoPseudofrequencies(),
    kargs...,
) where {T,A}
    sequences = getresiduesequences(msa)
    scores = sequencepairsmatrix(msa, T, Val{usediagonal}, diagonalvalue) # Named PairwiseListMatrix
    plm = getarray(scores)
    _mappairfreq!(
        f,
        sequences,
        plm,
        table,
        Val{usediagonal};
        weights = weights,
        pseudocounts = pseudocounts,
        pseudofrequencies = pseudofrequencies,
        kargs...,
    )
    scores
end

# DEPRECATED, usediagonal is now a boolean keyword argument
function mapseqpairfreq!(f, msa, table, usediagonal::Type{Val{D}}; kargs...) where {D}
    Base.depwarn(
        "The `usediagonal` positional argument of mapseqpairfreq! taking `Val{true}` or `Val{false}` is deprecated. Use `usediagonal = true` or `usediagonal = false` instead.",
        :mapseqpairfreq!,
        force = true,
    )
    mapseqpairfreq!(f, msa, table; usediagonal = D, kargs...)
end

# General mapfreq methods
# =======================

# The doctest for rank=2 can not be supported simultaneously for Julia 1.6 and 1.10
# because of the different behavior when printing the matrix headers.
# Therefore, the following docstest could be added once we drop support for Julia 1.6:
# mapfreq(sum, msa, rank=2)
# mapfreq(sum, msa, dims=1, rank=2)

"""
    mapfreq(f, msa; rank = 1, dims = 2, alphabet = UngappedAlphabet(), 
        weights = NoClustering(), pseudocounts = NoPseudocount(), 
        pseudofrequencies = NoPseudofrequencies(), probabilities = true, 
        usediagonal = false, diagonalvalue = NaN, kargs...)

It efficiently map a function (first argument) that takes a table of `Frequencies` or
`Probabilities` (depending on the `probabilities` keyword argument) calculated on 
sequences (`dims = 1`) or columns (`dims = 2`, the default) of an `msa` (second argument). 
If `rank = 1`, the default, the function is applied to each sequence or column. If 
`rank = 2`, the function is applied to each pair of sequences or columns. In that case, 
we can set the `usediagonal` keyword argument to `true` to apply the function to pairs
of the same sequence or column. The `diagonalvalue` keyword argument is used to set the 
value of the diagonal elements if `usediagonal` is `false`. By default, the function is not
applied to the diagonal elements (i.e. `usediagonal = false`) and the `diagonalvalue` is
set to `NaN`. The `alphabet` keyword argument can be used to set the alphabet used to 
construct the contingency table. The function also accepts the following keyword arguments:

$_mapfreq_kargs_doc

Note that the `pseudofrequencies` argument is only valid if `probabilities = true`. All the 
other keyword arguments are passed to the function `f`.
```jldoctest
julia> using Random, MIToS.MSA, MIToS.Information

julia> msa = rand(Random.MersenneTwister(1), Residue, 3, 6) # random MSA as an example
3×6 Matrix{Residue}:
 F  A  F  D  E  V
 T  R  R  G  F  I
 N  V  S  W  Q  T

julia> mapfreq(sum, msa) # default: rank=1, dims=2, probabilities=true
1×6 Named Matrix{Float64}
Function ╲ Col │   1    2    3    4    5    6
───────────────┼─────────────────────────────
sum            │ 1.0  1.0  1.0  1.0  1.0  1.0

julia> mapfreq(sum, msa, probabilities=false)
1×6 Named Matrix{Float64}
Function ╲ Col │   1    2    3    4    5    6
───────────────┼─────────────────────────────
sum            │ 3.0  3.0  3.0  3.0  3.0  3.0

julia> mapfreq(sum, msa, dims=1)
3×1 Named Matrix{Float64}
Seq ╲ Function │ sum
───────────────┼────
1              │ 1.0
2              │ 1.0
3              │ 1.0

```
"""
function mapfreq(
    f::Function,
    msa::AbstractArray{Residue};
    rank::Int = 1,
    dims::Int = 2,
    alphabet::ResidueAlphabet = UngappedAlphabet(),
    weights::WeightTypes = NoClustering(),
    pseudocounts::Pseudocount = NoPseudocount(),
    pseudofrequencies::Pseudofrequencies = NoPseudofrequencies(),
    probabilities::Bool = true,
    usediagonal::Bool = false,
    diagonalvalue::Float64 = NaN,
    kargs...,
)
    # Ensure that the keyword arguments are correct
    @assert dims == 1 || dims == 2 "The dimension must be 1 (sequences) or 2 (columns)."
    @assert rank == 1 || rank == 2 "The rank must be 1 (single sequences or columns) or 2 (pairs)."
    if pseudofrequencies !== NoPseudofrequencies()
        @assert probabilities "Set `probabilities = true` to use pseudofrequencies."
    end
    # Define the table to apply the function
    _table = ContingencyTable(Float64, Val{rank}, alphabet)
    table = probabilities ? Probabilities(_table) : Frequencies(_table)
    if rank == 1
        if dims == 1
            mapseqfreq!(
                f,
                msa,
                table;
                weights = weights,
                pseudocounts = pseudocounts,
                pseudofrequencies = pseudofrequencies,
                kargs...,
            )
        else
            mapcolfreq!(
                f,
                msa,
                table;
                weights = weights,
                pseudocounts = pseudocounts,
                pseudofrequencies = pseudofrequencies,
                kargs...,
            )
        end
    else
        if dims == 1
            mapseqpairfreq!(
                f,
                msa,
                table;
                usediagonal = usediagonal,
                diagonalvalue = diagonalvalue,
                weights = weights,
                pseudocounts = pseudocounts,
                pseudofrequencies = pseudofrequencies,
                kargs...,
            )
        else
            mapcolpairfreq!(
                f,
                msa,
                table;
                usediagonal = usediagonal,
                diagonalvalue = diagonalvalue,
                weights = weights,
                pseudocounts = pseudocounts,
                pseudofrequencies = pseudofrequencies,
                kargs...,
            )
        end
    end
end

# cMI
# ===

"""
`cumulative` allows to calculate cumulative scores (i.e. cMI) as defined
in  *Marino Buslje et al. 2010*:

> "We calculated a cumulative mutual information score (cMI) for each residue as the
> sum of MI values above a certain threshold for every amino acid pair where the particular
> residue appears. This value defines to what degree a given amino acid takes part in a
> mutual information network."

# References

  - [Marino Buslje, Cristina, et al. "Networks of high mutual information define the
    structural proximity of catalytic sites: implications for catalytic residue
    identification." PLoS computational biology 6.11 (2010):
    e1000978.](@cite 10.1371/journal.pcbi.1000978)
"""
function cumulative(plm::PairwiseListMatrix{T,D,VT}, threshold::T) where {T,D,VT}
    N = size(plm, 1)
    out = zeros(T, N)
    @iterateupper plm false begin
        elem = list[k]
        if !isnan(elem) && elem >= threshold
            out[i] += elem
            out[j] += elem
        end
    end
    reshape(out, (1, N))
end

function cumulative(
    nplm::NamedArray{T,2,PairwiseListMatrix{T,D,TV},DN},
    threshold::T,
) where {T,D,TV,DN}
    N = size(nplm, 2)
    name_list = names(nplm, 2)
    NamedArray(
        cumulative(getarray(nplm), threshold),
        (
            OrderedDict{String,Int}("cumulative" => 1),
            OrderedDict{String,Int}(name_list[i] => i for i = 1:N),
        ),
        ("Function", dimnames(nplm, 2)),
    )
end
