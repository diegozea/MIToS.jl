# Initialize an empty PairwiseListMatrix (PLM) for a pairwise measure in column or
# sequence pairs. These functions take care of the labels.
"""
Initialize an empty `PairwiseListMatrix` for a pairwise measure in column pairs. It uses the
column mapping (column number in the input MSA file) if it’s available, otherwise it uses
the actual column numbers. You can use the positional argument to indicate the number
`Type` (default: `Float64`), if the `PairwiseListMatrix` should store the diagonal values
on the list (default: `false`) and a default value for the diagonal (default: `NaN`).
"""
function sequencepairsmatrix{T}(msa::AbstractMatrix{Residue}, ::Type{T},
                                diagonal::Bool, diagonalvalue::T)
    plm = PairwiseListMatrix(T, nsequences(msa), diagonal, diagonalvalue)
    setlabels(plm, sequencenames(msa))
end

function sequencepairsmatrix(msa::AbstractMatrix{Residue})
    sequencepairsmatrix(msa, Float64, false, NaN)
end

"""
Initialize an empty `PairwiseListMatrix` for a pairwise measure in sequence pairs. It uses
the sequence names if they are available, otherwise it uses the actual sequence numbers.
You can use the positional argument to indicate the number `Type` (default: `Float64`),
if the `PairwiseListMatrix` should store the diagonal values on the list (default: `false`)
and a default value for the diagonal (default: `NaN`).
"""
function columnpairsmatrix{T}(msa::AbstractMatrix{Residue}, ::Type{T},
                              diagonal::Bool, diagonalvalue::T)
    plm = PairwiseListMatrix(T, ncolumns(msa), diagonal, diagonalvalue)
    setlabels(plm, columnnames(msa))
end

function columnpairsmatrix(msa::AbstractMatrix{Residue})
    columnpairsmatrix(msa, Float64, false, NaN)
end
