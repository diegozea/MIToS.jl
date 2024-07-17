for T in (
    :(AlignedSequence),
    :(AnnotatedSequence),
    :(AnnotatedAlignedSequence),
    :(MultipleSequenceAlignment),
    :(AnnotatedMultipleSequenceAlignment),
)
    @eval Base.IndexStyle(::Type{$(T)}) = Base.IndexLinear()
end

# ---
# Performant lookup of the sequence and column positions.
# CAUTION: This relies on the NamedArrays internal representation.

# This is a bit hacky, but it's the fastest way to get the index.
"""
    sequence_index(msa, seq_name)

Return the index (integer position) of the sequence with name `seq_name` in the MSA `msa`.
A `KeyError` is thrown if the sequence name does not exist. If `seq_name` is an integer,
the same integer is returned without checking if it is a valid index.
"""
sequence_index(msa::NamedResidueMatrix, seq_name::AbstractString) = msa.dicts[1][seq_name]

"""
    column_index(msa, col_name)

Return the index (integer position) of the column with name `col_name` in the MSA `msa`.
A `KeyError` is thrown if the column name does not exist. If `col_name` is an integer,
the same integer is returned without checking if it is a valid index.
"""
column_index(msa::NamedResidueMatrix, col_name::AbstractString) = msa.dicts[2][col_name]

function sequence_index(msa::AbstractResidueMatrix, seq_name::AbstractString)
    sequence_index(msa.matrix, seq_name)
end

function column_index(msa::AbstractResidueMatrix, column_name::AbstractString)
    column_index(msa.matrix, column_name)
end

# Do not allow indexing Matrix{Residue} with AbstractString

function sequence_index(msa::Matrix{Residue}, seq_name::AbstractString)
    throw(
        ErrorException(
            "There are no sequence names in a Matrix{Residue} object, use an integer index instead.",
        ),
    )
end

function column_index(msa::Matrix{Residue}, column_name::AbstractString)
    throw(
        ErrorException(
            "There are no column names in a Matrix{Residue} object, use an integer index instead.",
        ),
    )
end

# If the user already gives a position, return the same position.

sequence_index(msa, sequence_index::Int) = sequence_index
column_index(msa, column_index::Int) = column_index

# ---

@inline Base.getindex(x::AbstractResidueMatrix, args...) = getindex(namedmatrix(x), args...)

@inline function Base.setindex!(x::AbstractResidueMatrix, value, args...)
    setindex!(namedmatrix(x), value, args...)
end

# Special getindex/setindex! for sequences to avoid `seq["seqname","colname"]`

@inline function Base.getindex(
    x::Union{AbstractSequence,AbstractAlignedSequence},
    i::Union{Int,AbstractString},
)
    getindex(namedmatrix(x), 1, i)
end

@inline function Base.setindex!(
    x::Union{AbstractSequence,AbstractAlignedSequence},
    value,
    i::Union{Int,AbstractString},
)
    setindex!(namedmatrix(x), value, 1, i)
end

# ## Special getindex to conserve annotations and types

"""
Helper function to create a boolean mask to select sequences.
"""
function _get_selected_sequences(msa, selector)
    type = eltype(selector)
    if type == Bool
        return selector
    else
        to_select = Set(selector)
        if type <: Number
            Bool[i in to_select for i = 1:nsequences(msa)]
        elseif type <: AbstractString
            Bool[i in to_select for i in sequencenames(msa)]
        else
            throw(ArgumentError("$type is not a valid element type for the selector."))
        end
    end
end

function _column_indices(msa, selector)
    if eltype(selector) <: AbstractString
        return Int[column_index(msa, i) for i in selector]
    end
    selector
end

function Base.getindex(
    msa::AnnotatedMultipleSequenceAlignment,
    seqs::AbstractArray,
    cols::AbstractArray,
)
    annot = copy(annotations(msa))
    seq_selector = _get_selected_sequences(msa, seqs)
    filtersequences!(annot, sequencenames(msa), seq_selector)
    _annotate_seq_modification!(annot, seq_selector)
    col_selector = _column_indices(msa, cols)
    filtercolumns!(annot, col_selector)
    _annotate_col_modification!(annot, col_selector)
    AnnotatedMultipleSequenceAlignment(msa.matrix[seqs, cols], annot)
end

function Base.getindex(
    msa::AnnotatedMultipleSequenceAlignment,
    seqs::AbstractArray,
    cols::Colon,
)
    annot = copy(annotations(msa))
    seq_selector = _get_selected_sequences(msa, seqs)
    filtersequences!(annot, sequencenames(msa), seq_selector)
    _annotate_seq_modification!(annot, seq_selector)
    AnnotatedMultipleSequenceAlignment(msa.matrix[seqs, cols], annot)
end

function Base.getindex(
    msa::AnnotatedMultipleSequenceAlignment,
    seqs::Colon,
    cols::AbstractArray,
)
    annot = copy(annotations(msa))
    col_selector = _column_indices(msa, cols)
    filtercolumns!(annot, col_selector)
    _annotate_col_modification!(annot, col_selector)
    AnnotatedMultipleSequenceAlignment(msa.matrix[seqs, cols], annot)
end

Base.getindex(msa::AnnotatedMultipleSequenceAlignment, seqs::Colon, cols::Colon) = copy(msa)

function Base.getindex(
    msa::MultipleSequenceAlignment,
    seqs::Union{AbstractArray,Colon},
    cols::Union{AbstractArray,Colon},
)
    MultipleSequenceAlignment(msa.matrix[seqs, cols])
end

function Base.getindex(seq::AnnotatedAlignedSequence, cols::AbstractArray)
    seq_copy = copy(seq)
    col_selector = _column_indices(seq, cols)
    filtercolumns!(seq_copy, col_selector)
    _annotate_col_modification!(seq_copy, col_selector)
    seq_copy
end

# TODO: AnnotatedSequence

function Base.getindex(seq::AlignedSequence, cols::Union{AbstractArray,Colon})
    AlignedSequence(seq.matrix[cols])
end
