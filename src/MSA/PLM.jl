# Initialize an empty PairwiseListMatrix (PLM) for a pairwise measure in column or sequence pairs.
# These functions take care of the labels.

"Returns the column labels: use the column mapping (column number in the input MSA file) if it’s available, otherwise it uses the actual column numbers."
columnlabels(msa::AnnotatedMultipleSequenceAlignment) = haskey(getannotfile(msa), "ColMap") ? getcolumnmapping(msa) : collect(1:ncolumns(msa))
columnlabels(msa::AbstractMatrix{Residue}) = collect(1:ncolumns(msa))

"Returns the sequence labels: use the sequence names if they are available, otherwise it uses the actual sequence numbers."
sequencelabels(msa::Union{AnnotatedMultipleSequenceAlignment, MultipleSequenceAlignment}) = ASCIIString[ name for name in names(msa) ]
sequencelabels(msa::AbstractMatrix{Residue}) = ASCIIString[ string(number) for number in collect(1:nsequences(msa)) ]

"""
Initialize an empty `PairwiseListMatrix` for a pairwise measure in column pairs.
This function uses the column mapping (column number in the input MSA file)
if it’s available, otherwise it uses the actual column numbers.
You can use the positional argument to indicate the number `Type` (default: `Float64`),
if the `PairwiseListMatrix` should store the diagonal values on the list (default: `false`) and
a default value for the diagonal (default: `NaN`).
"""
function sequencepairsmatrix{T}(msa::AbstractMatrix{Residue}, ::Type{T}, diagonal::Bool, diagonalvalue::T)
    PairwiseListMatrix(T, nsequences(msa), sequencelabels(msa), diagonal, diagonalvalue)
end

sequencepairsmatrix(msa::AbstractMatrix{Residue}) = sequencepairsmatrix(msa, Float64, false, NaN)

"""
Initialize an empty `PairwiseListMatrix` for a pairwise measure in sequence pairs.
This function uses the sequence names if they are available, otherwise it uses the actual sequence numbers.
You can use the positional argument to indicate the number `Type` (default: `Float64`),
if the `PairwiseListMatrix` should store the diagonal values on the list (default: `false`) and
a default value for the diagonal (default: `NaN`).
"""
function columnpairsmatrix{T}(msa::AbstractMatrix{Residue}, ::Type{T}, diagonal::Bool, diagonalvalue::T)
    PairwiseListMatrix(T, ncolumns(msa), columnlabels(msa), diagonal, diagonalvalue)
end

columnpairsmatrix(msa::AbstractMatrix{Residue}) = columnpairsmatrix(msa, Float64, false, NaN)
