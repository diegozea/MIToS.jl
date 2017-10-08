# Initialize an empty PairwiseListMatrix (PLM) for a pairwise measure in column or
# sequence pairs. These functions take care of the labels.
"""
Initialize an empty `PairwiseListMatrix` for a pairwise measure in column pairs. It uses the
column mapping (column number in the input MSA file) if it’s available, otherwise it uses
the actual column numbers. You can use the positional argument to indicate the number
`Type` (default: `Float64`), if the `PairwiseListMatrix` should store the diagonal values
on the list (default: `false`) and a default value for the diagonal (default: `NaN`).
"""
function sequencepairsmatrix(msa::AbstractMatrix{Residue}, ::Type{T},
                             ::Type{Val{diagonal}}, diagonalvalue::T) where {T,diagonal}
    plm = PairwiseListMatrix(T, nsequences(msa), diagonal, diagonalvalue)
    nplm = setlabels(plm, sequencenames(msa))
    setdimnames!(nplm, ["Seq1", "Seq2"])
    nplm::NamedArray{T,2,PairwiseListMatrix{T,diagonal,Vector{T}},NTuple{2,OrderedDict{String,Int}}}
end

function sequencepairsmatrix(msa::AbstractMatrix{Residue})
    sequencepairsmatrix(msa, Float64, Val{false}, NaN)
end

"""
Initialize an empty `PairwiseListMatrix` for a pairwise measure in sequence pairs. It uses
the sequence names if they are available, otherwise it uses the actual sequence numbers.
You can use the positional argument to indicate the number `Type` (default: `Float64`),
if the `PairwiseListMatrix` should store the diagonal values on the list (default: `false`)
and a default value for the diagonal (default: `NaN`).
"""
function columnpairsmatrix(msa::AbstractMatrix{Residue}, ::Type{T},
                           ::Type{Val{diagonal}}, diagonalvalue::T) where {T,diagonal}
    plm = PairwiseListMatrix(T, ncolumns(msa), diagonal, diagonalvalue)
    nplm = setlabels(plm, columnnames(msa))
    setdimnames!(nplm, ["Col1", "Col2"])
    nplm::NamedArray{T,2,PairwiseListMatrix{T,diagonal,Vector{T}},NTuple{2,OrderedDict{String,Int}}}
end

function columnpairsmatrix(msa::AbstractMatrix{Residue})
    columnpairsmatrix(msa, Float64, Val{false}, NaN)
end
