for T in (:(AlignedSequence),
    :(AnnotatedAlignedSequence),
    :(MultipleSequenceAlignment),
    :(AnnotatedMultipleSequenceAlignment))
    @eval Base.IndexStyle(::Type{$(T)}) = Base.IndexLinear()
end

@inline Base.getindex(x::AbstractAlignedObject,
              args...) = getindex(namedmatrix(x), args...)

@inline function Base.setindex!(x::AbstractAlignedObject, value, args...)
    setindex!(namedmatrix(x), value, args...)
end

# Special getindex/setindex! for sequences to avoid `seq["seqname","colname"]`

@inline function Base.getindex(x::AbstractAlignedSequence,
                       i::Union{Int,AbstractString})
    getindex(namedmatrix(x), 1, i)
end

@inline function Base.setindex!(x::AbstractAlignedSequence, value,
                        i::Union{Int,AbstractString})
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
        if type  <: Number
            Bool[i in to_select for i in 1:nsequences(msa)]
        elseif type <: AbstractString
            Bool[i in to_select for i in sequencenames(msa)]
        end
    end
end

_get_names_to_index(msa::AbstractAlignedObject) = NamedArrays.names(msa.matrix, 2)
_get_names_to_index(msa::NamedResidueMatrix{T}) where T = NamedArrays.names(msa, 2)
function _get_names_to_index(msa::Matrix{Residue})
    throw(ErrorException("`Matrix{Residue}` cannot be indexed with `AbstractString`s."))
end

function _column_indices(msa, selector)
    if eltype(selector) <: AbstractString
        colnames = _get_names_to_index(msa)
        return Int[findfirst(==(i), colnames) for i in selector]
    end
    selector
end

function Base.getindex(msa::AnnotatedMultipleSequenceAlignment, 
               seqs::AbstractArray, cols::AbstractArray)
    annot = copy(annotations(msa))
    seq_selector = _get_selected_sequences(msa, seqs)
    filtersequences!(annot, sequencenames(msa), seq_selector)
    _annotate_seq_modification!(annot, seq_selector)
    col_selector = _column_indices(msa, cols)
    filtercolumns!(annot, col_selector)
    _annotate_col_modification!(annot, col_selector)
    AnnotatedMultipleSequenceAlignment(msa.matrix[seqs, cols], annot)
end

function Base.getindex(msa::AnnotatedMultipleSequenceAlignment, 
               seqs::AbstractArray, cols::Colon)
    annot = copy(annotations(msa))
    seq_selector = _get_selected_sequences(msa, seqs)
    filtersequences!(annot, sequencenames(msa), seq_selector)
    _annotate_seq_modification!(annot, seq_selector)
    AnnotatedMultipleSequenceAlignment(msa.matrix[seqs, cols], annot)
end

function Base.getindex(msa::AnnotatedMultipleSequenceAlignment, 
               seqs::Colon, cols::AbstractArray)
    annot = copy(annotations(msa))
    col_selector = _column_indices(msa, cols)
    filtercolumns!(annot, col_selector)
    _annotate_col_modification!(annot, col_selector)
    AnnotatedMultipleSequenceAlignment(msa.matrix[seqs, cols], annot)
end

Base.getindex(msa::AnnotatedMultipleSequenceAlignment, seqs::Colon, cols::Colon) = copy(msa)

function Base.getindex(msa::MultipleSequenceAlignment, 
seqs::Union{AbstractArray,Colon}, 
cols::Union{AbstractArray,Colon})
    MultipleSequenceAlignment(msa.matrix[seqs, cols])
end

function Base.getindex(seq::AnnotatedAlignedSequence, cols::AbstractArray)
    seq_copy = copy(seq)
    col_selector = _column_indices(msa, cols)
    filtercolumns!(seq_copy, col_selector)
    _annotate_col_modification!(seq_copy, col_selector)
    seq_copy
end

function Base.getindex(seq::AlignedSequence, cols::Union{AbstractArray,Colon})
    AlignedSequence(seq.matrix[cols])
end
