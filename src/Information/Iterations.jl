_get_matrix_residue(msa::AbstractMultipleSequenceAlignment) = NamedArrays.array(namedmatrix(msa))
_get_matrix_residue(msa::NamedArray) = NamedArrays.array(msa)
_get_matrix_residue(msa) = msa

# Kernel function to fill ContingencyTable based on residues
function _mapfreq_kernel!{T,N,A}(f,
                                table::ContingencyTable{T,N,A},
                                probabilities,
                                weights, pseudocounts, pseudofrequencies,
                                res::NTuple{N,AbstractVector{Residue}})
    _cleanup_table!(table) # count! call _cleanup_temporal! and cleans marginals
    count!(table, weights, pseudocounts, res...)
    apply_pseudocount!(table, pseudocounts)
    _update!(table)
    if probabilities
        normalize!(table)
        apply_pseudofrequencies!(table, pseudofrequencies)
    end
    f(table)
end

_mapfreq_kargs_doc = """
- `probabilities` (default: `true`): Indicates if the table should be normalized before applying the function.
- `weights` (default: `NoClustering()`): Weights to be used for table counting.
- `pseudocounts` (default: `NoPseudocount()`): `Pseudocount` object to be applied to table.
- `pseudofrequencies` (default: `NoPseudofrequencies()`): If `probabilities` is `true`, `Pseudofrequencies` to be applied to the normalized table.
"""

_mappairfreq_kargs_doc = """
- `usediagonal` (default: `true`): Indicates if the function should be applied to identical element pairs.
- `diagonalvalue` (default: `0`): Value to fill diagonal elements if `usediagonal` is `false`.
"""

# Residues: The output is a Named Vector
# --------------------------------------

function _mapfreq!{T,A,V<:AbstractArray{Residue}}(f::Function,
                        res_list::Vector{V}, # sequences or columns
                        table::ContingencyTable{T,1,A};
                        probabilities::Bool = true,
                        weights = NoClustering(),
                        pseudocounts::Pseudocount = NoPseudocount(),
                        pseudofrequencies::Pseudofrequencies = NoPseudofrequencies())
    scores = map(res_list) do res
        _mapfreq_kernel!(f, table, probabilities,
                         weights, pseudocounts, pseudofrequencies, (res,))
    end
    scores
end

# Map to each column

"""
$_mapfreq_kargs_doc
"""
function mapcolfreq!{T,A}(f::Function, msa::AbstractMatrix{Residue},
                          table::ContingencyTable{T,1,A}; kargs...)
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
$_mapfreq_kargs_doc
"""
function mapseqfreq!{T,A}(f::Function, msa::AbstractMatrix{Residue},
                          table::ContingencyTable{T,1,A}; kargs...)
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

function _mappairfreq!{T,D,TV,A,V<:AbstractArray{Residue}}(f::Function,
                                res_list::Vector{V}, # sequences or columns
                                plm::PairwiseListMatrix{T,D,TV}, # output
                                table::ContingencyTable{T,2,A},
                                usediagonal::Type{Val{D}} = Val{true};
                                probabilities::Bool = true,
                                weights = NoClustering(),
                                pseudocounts::Pseudocount = NoPseudocount(),
                                pseudofrequencies::Pseudofrequencies = NoPseudofrequencies(),
                                diagonalvalue::T = zero(T))
    @inbounds @iterateupper plm D begin
        list[k] = :($_mapfreq_kernel!)(:($f), :($table), :($probabilities),
                                       :($weights),
                                       :($pseudocounts), :($pseudofrequencies),
                                       (:($res_list)[i], :($res_list)[j]))
    end
    plm
end

# Map to column pairs

"""
$_mapfreq_kargs_doc
$_mappairfreq_kargs_doc
"""
function mapcolpairfreq!{T,A,D}(f::Function, msa::AbstractMatrix{Residue},
                                table::ContingencyTable{T,2,A},
                                ::Type{Val{D}} = Val{true};
                                diagonalvalue::T = zero(T),
                                kargs...)
    ncol = ncolumns(msa)
    residues = _get_matrix_residue(msa)
    columns = map(i -> view(residues,:,i), 1:ncol) # 2x faster than calling view inside the loop
    scores = columnpairsmatrix(msa, T, Val{D}, diagonalvalue) # Named PairwiseListMatrix
    plm = NamedArrays.array(scores) # PairwiseListMatrix
    _mappairfreq!(f, columns, plm, table, Val{D};
                  diagonalvalue=diagonalvalue, kargs...)
    scores
end

# Map to sequence pairs

"""
$_mapfreq_kargs_doc
$_mappairfreq_kargs_doc
"""
function mapseqpairfreq!{T,A,D}(f::Function, msa::AbstractMatrix{Residue},
                                table::ContingencyTable{T,2,A},
                                ::Type{Val{D}} = Val{true};
                                diagonalvalue::T = zero(T),
                                kargs...)
    sequences = getresiduesequences(msa)
    scores = sequencepairsmatrix(msa, T, Val{D}, diagonalvalue) # Named PairwiseListMatrix
    plm = NamedArrays.array(scores) # PairwiseListMatrix
    _mappairfreq!(f, sequences, plm, table, Val{D};
                  diagonalvalue=diagonalvalue, kargs...)
    scores
end

# TO DO:
# N >= 3 using using combinations from Combinatorics.jl

# """
# `Information.mapcolfreq!(aln, [count,] use, [α, β,] measure, [pseudocount,] [weight,] [usediagonal, diagonalvalue])`
#
# This function `estimate` a `AbstractMeasure` in columns or pair of columns of a MSA.
#
# - `aln` : This argument is mandatory and it should be a `Matrix{Residue}`. Use the function `getresidues` (from the MSA module) over a MSA object to get the needed matrix.
# - `count` : It should be defined when `use` is a `ResidueProbability` object. It indicates the element type of the counting table.
# - `use` : This argument is mandatory and indicates the sub-type of `ResidueContingencyTables` used by `estimate` inside the function.
# If the table has one dimension (`N`=`1`), the occurrences/probabilities are counted for each sequence/column.
# If the table has two dimension (`N`=`2`), pairs of sequences/columns are used.
# The dimension `N` and the `UseGap` parameter of `Residueount{T, N, UseGap}` or `ResidueProbability{T, N, UseGap}` determines the output and behaviour of this functions.
# If `UseGap` is true, gaps are used in the estimations.
# - `α` : This argument is optional, and indicates the weight of real frequencies to apply BLOSUM62 based pseudo frequencies.
# - `β` : This argument is optional, and indicates the weight of BLOSUM62 based pseudo frequencies.
# - `measure` : This argument is mandatory and indicates the measure to be used by `estimate` inside the function.
# - `pseudocount` : This argument is optional. It should be an `AdditiveSmoothing` instance (default to zero).
# - `weight` : This argument is optional. It should be an instance of `ClusteringResult` or `AbstractVector` (vector of weights).
# Each sequence has weight 1 (`NoClustering()`) by default.
# - `usediagonal` : This functions return a `Vector` in the one dimensional case, or a `PairwiseListMatrix` in the bidimensional case.
# This argument only have sense in the bidimensional case and indicates if the list on the `PairwiseListMatrix` should include the diagonal (default to `true`).
# - `diagonalvalue` : This argument is optional (default to zero). Indicates the value of output diagonal elements.
# """

# cMI
# ===

"""
`cumulative` allows to calculate cumulative scores (i.e. cMI) as defined in Buslje et. al. 2010

*"We calculated a cumulative mutual information score (cMI) for each residue as the sum of
MI values above a certain threshold for every amino acid pair where the particular residue
appears. This value defines to what degree a given amino acid takes part in a mutual
information network."*
Buslje, Cristina Marino, Elin Teppa, Tomas Di Doménico, José María Delfino, and Morten Nielsen.
*Networks of high mutual information define the structural proximity of catalytic sites:
implications for catalytic residue identification.*
PLoS Comput Biol 6, no. 11 (2010): e1000978.
"""
function cumulative{T,D,VT}(plm::PairwiseListMatrix{T,D,VT}, threshold::T)
    N = size(nplm, 1)
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

function cumulative{T,D,TV,DN}(nplm::NamedArray{T,2,PairwiseListMatrix{T,D,TV},DN},
                               threshold::T)
    N = size(nplm, 2)
    name_list = names(nplm, 2)
    NamedArray(cumulative(NamedArrays.array(nplm), threshold),
               (OrderedDict{String,Int}("cumulative" => 1),
               OrderedDict{String,Int}(name_list[i] => i for i in 1:N)),
               ("Function", dimnames(nplm,2)))
end
