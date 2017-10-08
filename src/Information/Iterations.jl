function _get_matrix_residue(msa::AbstractMultipleSequenceAlignment)::Matrix{Residue}
    getarray(namedmatrix(msa))
end
_get_matrix_residue(msa::NamedArray)::Matrix{Residue} = getarray(msa)
_get_matrix_residue(msa::Matrix{Residue})::Matrix{Residue} = msa

# Kernel function to fill ContingencyTable based on residues
function _mapfreq_kernel!(f::Function,
                          freqtable::Union{Probabilities{T,N,A},Counts{T,N,A}},
                          weights, pseudocounts, pseudofrequencies,
                          res::Vararg{AbstractVector{Residue},N}) where {T,N,A}
    table = freqtable.table
    _cleanup_table!(table) # count! calls _cleanup_temporal! and cleans marginals
    count!(table, weights, pseudocounts, res...) # count! calls apply_pseudocount! and  _update!
    if isa(freqtable, Probabilities{T,N,A})
        normalize!(table)
        apply_pseudofrequencies!(table, pseudofrequencies)
    end
    f(freqtable)
end

_mapfreq_kargs_doc = """
- `weights` (default: `NoClustering()`): Weights to be used for table counting.
- `pseudocounts` (default: `NoPseudocount()`): `Pseudocount` object to be applied to table.
- `pseudofrequencies` (default: `NoPseudofrequencies()`): `Pseudofrequencies` to be applied to the normalized (probabilities) table.
"""

_mappairfreq_kargs_doc = """
- `diagonalvalue` (default: `0`): Value to fill diagonal elements if `usediagonal` is `Val{false}`.
"""

# Residues: The output is a Named Vector
# --------------------------------------

function _mapfreq!(f::Function,
res_list::Vector{V}, # sequences or columns
table::Union{Probabilities{T,1,A},Counts{T,1,A}};
weights = NoClustering(),
pseudocounts::Pseudocount = NoPseudocount(),
pseudofrequencies::Pseudofrequencies = NoPseudofrequencies()) where {T,A,V<:AbstractArray{Residue}}
    scores = map(res_list) do res
        _mapfreq_kernel!(f, table, weights, pseudocounts, pseudofrequencies, res)
    end
    scores
end

# Map to each column

"""
It efficiently map a function (first argument) that takes a table of `Counts` or
`Probabilities` (third argument). The table is filled in place with the counts or
probabilities of each column from the `msa` (second argument).

$_mapfreq_kargs_doc
"""
function mapcolfreq!(f::Function, msa::AbstractMatrix{Residue},
                     table::Union{Probabilities{T,1,A},Counts{T,1,A}}; kargs...) where {T,A}
    N = ncolumns(msa)
    residues = _get_matrix_residue(msa)
    ncol = size(residues,2)
    columns = map(i -> view(residues,:,i), 1:ncol) # 2x faster than calling view inside the loop
    scores = _mapfreq!(f, columns, table; kargs...)
    name_list = columnnames(msa)
    NamedArray(reshape(scores, (1,N)),
               (OrderedDict{String,Int}(string(f) => 1),
                OrderedDict{String,Int}(name_list[i] => i for i in 1:N)),
               ("Function","Col"))
end

# Map to each sequence

"""
It efficiently map a function (first argument) that takes a table of `Counts` or
`Probabilities` (third argument). The table is filled in place with the counts or
probabilities of each sequence from the `msa` (second argument).

$_mapfreq_kargs_doc
"""
function mapseqfreq!(f::Function, msa::AbstractMatrix{Residue},
                     table::Union{Probabilities{T,1,A},Counts{T,1,A}}; kargs...) where {T,A}
    N = nsequences(msa)
    sequences = getresiduesequences(msa)
    name_list = sequencenames(msa)
    scores = _mapfreq!(f, sequences, table; kargs...)
    NamedArray(reshape(scores, (N,1)),
               (OrderedDict{String,Int}(name_list[i] => i for i in 1:N),
               OrderedDict{String,Int}(string(f) => 1)),
               ("Seq","Function"))
end

# Residue pairs: The output is a Named PairwiseListMatrix
# -------------------------------------------------------

function _mappairfreq!(f::Function,
res_list::Vector{V}, # sequences or columns
plm::PairwiseListMatrix{T,D,TV}, # output
table::Union{Probabilities{T,2,A},Counts{T,2,A}},
usediagonal::Type{Val{D}} = Val{true};
weights = NoClustering(),
pseudocounts::Pseudocount = NoPseudocount(),
pseudofrequencies::Pseudofrequencies = NoPseudofrequencies(),
diagonalvalue::T = zero(T)) where {T,D,TV,A,V<:AbstractArray{Residue}} # diagonalvalue used by mapcolpairfreq! ...
    @inbounds @iterateupper plm D begin
        list[k] = :($_mapfreq_kernel!)(:($f), :($table), :($weights),
                                       :($pseudocounts), :($pseudofrequencies),
                                       :($res_list)[i], :($res_list)[j])
    end
    plm
end

# Map to column pairs

"""
It efficiently map a function (first argument) that takes a table of `Counts` or
`Probabilities` (third argument). The table is filled in place with the counts or
probabilities of each pair of columns from the `msa` (second argument).
The fourth positional argument `usediagonal` indicates if the function should be applied
to identical element pairs (default to `Val{true}`).

$_mapfreq_kargs_doc
$_mappairfreq_kargs_doc
"""
function mapcolpairfreq!(f::Function, msa::AbstractMatrix{Residue},
                         table::Union{Probabilities{T,2,A},Counts{T,2,A}},
                         usediagonal::Type{Val{D}} = Val{true};
                         diagonalvalue::T = zero(T),
                         kargs...) where {T,A,D}
    ncol = ncolumns(msa)
    residues = _get_matrix_residue(msa)
    columns = map(i -> view(residues,:,i), 1:ncol) # 2x faster than calling view inside the loop
    scores = columnpairsmatrix(msa, T, Val{D}, diagonalvalue) # Named PairwiseListMatrix
    plm = getarray(scores)::PairwiseListMatrix{T,D,Vector{T}} # PairwiseListMatrix
    _mappairfreq!(f, columns, plm, table, Val{D}; kargs...)
    scores
end

# Map to sequence pairs

"""
It efficiently map a function (first argument) that takes a table of `Counts` or
`Probabilities` (third argument). The table is filled in place with the counts or
probabilities of each pair of sequences from the `msa` (second argument).
The fourth positional argument `usediagonal` indicates if the function should be applied
to identical element pairs (default to `Val{true}`).

$_mapfreq_kargs_doc
$_mappairfreq_kargs_doc
"""
function mapseqpairfreq!(f::Function, msa::AbstractMatrix{Residue},
                         table::Union{Probabilities{T,2,A},Counts{T,2,A}},
                         usediagonal::Type{Val{D}} = Val{true};
                         diagonalvalue::T = zero(T),
                         kargs...) where {T,A,D}
    sequences = getresiduesequences(msa)
    scores = sequencepairsmatrix(msa, T, Val{D}, diagonalvalue) # Named PairwiseListMatrix
    plm = getarray(scores)::PairwiseListMatrix{T,D,Vector{T}} # PairwiseListMatrix
    _mappairfreq!(f, sequences, plm, table, Val{D}; kargs...)
    scores
end

# cMI
# ===

"""
`cumulative` allows to calculate cumulative scores (i.e. cMI) as defined in Buslje et. al. 2010

*"We calculated a cumulative mutual information score (cMI) for each residue as the sum of
MI values above a certain threshold for every amino acid pair where the particular residue
appears. This value defines to what degree a given amino acid takes part in a mutual
information network."*
Buslje, Cristina Marino, Elin Teppa, Tomas Di Doménico, José María Delfino, and Morten
Nielsen. *Networks of high mutual information define the structural proximity of catalytic
sites: implications for catalytic residue identification.* PLoS Comput Biol 6, no. 11
(2010): e1000978.
"""
function cumulative(plm::PairwiseListMatrix{T,D,VT}, threshold::T) where {T,D,VT}
    N = size(plm, 1)
    out = zeros(T, N)
    @iterateupper plm false begin
        elem = list[k]
        if !isnan(elem) && elem >= :($threshold)
            :($out)[i] += elem
            :($out)[j] += elem
        end
    end
    reshape(out, (1,N))
end

function cumulative(nplm::NamedArray{T,2,PairwiseListMatrix{T,D,TV},DN},
                    threshold::T) where {T,D,TV,DN}
    N = size(nplm, 2)
    name_list = names(nplm, 2)
    NamedArray(cumulative(getarray(nplm), threshold),
               (OrderedDict{String,Int}("cumulative" => 1),
               OrderedDict{String,Int}(name_list[i] => i for i in 1:N)),
               ("Function", dimnames(nplm,2)))
end
