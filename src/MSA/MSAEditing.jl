# Filters
# =======

# Creation of a boolean mask array
# --------------------------------

# If the input is a matrix, dropdims the singleton dimension

function _column_mask(mask::AbstractMatrix{Bool}, msa)
    @assert size(mask, 1) == 1 "The mask should be a vector or a matrix of size (1,ncol)"
    @assert size(mask, 2) == ncolumns(msa) "One boolean value per column is needed."
    dropdims(mask, dims=1)
end

function _sequence_mask(mask::AbstractMatrix{Bool}, msa)
    @assert size(mask, 2) == 1 "The mask should be a vector or a matrix of size (nseq,1)"
    @assert size(mask, 1) == nsequences(msa) "One boolean value per sequence is needed."
    dropdims(mask, dims=2)
end

# If the input is a function, apply it to the respective slice

function _column_mask(mask::Function, msa)
    out = mapslices(mask, msa, dims=1)
    if size(out) != (1, ncolumns(msa)) || eltype(out) != Bool
        error("The function must return a Bool element per column.")
    end
    dropdims(out, dims=1)
end

function _sequence_mask(mask::Function, msa)
    out = mapslices(mask, msa, dims=2)
    if size(out) != (nsequences(msa), 1) || eltype(out) != Bool
        error("The function must return a Bool element per column.")
    end
    dropdims(out, dims=2)
end

# If the mask is an AbstractVector{Bool}, return it without changes

function _column_mask(mask::AbstractVector{Bool}, msa)
    @assert length(mask) == ncolumns(msa) "One boolean value per column is needed."
    mask
end

function _sequence_mask(mask::AbstractVector{Bool}, msa)
    @assert length(mask) == nsequences(msa) "One boolean value per sequence is needed."
    mask
end

# Filter sequences
# ----------------

"It's similar to `filtersequences!` but for an `AbstractMatrix{Residue}`"
filtersequences(msa::AbstractMatrix{Residue}, mask) = msa[_sequence_mask(mask, msa), :]

"""
`filtersequences!(msa, mask[, annotate::Bool=true])`

It allows to filter `msa` sequences using a `AbstractVector{Bool}` `mask`
(It removes sequnces with `false` values). `AnnotatedMultipleSequenceAlignment` annotations
are updated if `annotate` is `true` (default).
"""
function filtersequences!(msa::AnnotatedMultipleSequenceAlignment,
                          mask, annotate::Bool=true)
    boolean_vector = _sequence_mask(mask, msa)
    filtersequences!(annotations(msa), sequencenames(msa), boolean_vector)
    # Filter annotations first, names will be changed by filtersequences:
    msa.matrix = filtersequences(namedmatrix(msa), boolean_vector)
    annotate && annotate_modification!(msa, string("filtersequences! : ",
                                                    sum(.~boolean_vector),
                                                    " sequences have been deleted."))
    msa
end

function filtersequences!(msa::MultipleSequenceAlignment,
                          mask, annotate::Bool=false) # annotate is useful for calling this
                                                      # inside other functions
    msa.matrix = filtersequences(namedmatrix(msa), _sequence_mask(mask, msa))
    msa
end

function filtersequences(msa::Union{AnnotatedMultipleSequenceAlignment,
                                    MultipleSequenceAlignment}, args...)
    filtersequences!(deepcopy(msa), args...)
end

# Filter columns
# --------------

"It's similar to `filtercolumns!` but for an `AbstractMatrix{Residue}`"
function filtercolumns(msa::AbstractMatrix{Residue}, mask)
    msa[:, _column_mask(mask, msa)]
end

"""
`filtercolumns!(msa, mask[, annotate::Bool=true])`

It allows to filter MSA or aligned sequence columns/positions using a
`AbstractVector{Bool}` `mask`. Annotations are updated if `annotate` is `true` (default).
"""
function filtercolumns!(x::AnnotatedAlignedObject,
                        mask, annotate::Bool=true)
    boolean_vector = _column_mask(mask, x)
    filtercolumns!(annotations(x), boolean_vector)
    x.matrix = filtercolumns(namedmatrix(x), boolean_vector)
    annotate && annotate_modification!(x,string("filtercolumns! : ", sum(.~boolean_vector),
                                                " columns have been deleted."))
    x
end

function filtercolumns!(x::UnannotatedAlignedObject,
                        mask, annotate::Bool=false) # annotate is useful for calling this
                                                    # inside other functions
    x.matrix = filtercolumns(namedmatrix(x), _column_mask(mask, x))
    x
end

filtercolumns(x::AbstractAlignedObject, args...) = filtercolumns!(deepcopy(x), args...)

# Util function
# -------------

"""
This function takes a vector of sequence names and a sequence id.
It returns the position of that id in the vector.
If the id isn't in the vector, It throws an error.
"""
function _get_seqid_index(names::Vector{String}, sequence_id::String)
    id_index = something(findfirst(isequal(sequence_id), names), 0)
    if id_index == 0
        throw(ErrorException("$sequence_id is not in the list of sequence names."))
    end
    id_index
end

# Reference
# ---------
"It swaps the names on the positions `i` and `j` of a `Vector{String}`"
function _swap!(names::Vector{String}, i::Int, j::Int)
    name = names[i,:]
    names[i,:] = names[j,:]
    names[j,:] = name
    names
end

"""
It swaps the sequences on the positions `i` and `j` of an MSA. Also it's possible to swap
sequences using their sequence names/identifiers when the MSA object as names.
"""
function swapsequences!(matrix::Matrix{Residue}, i::Int, j::Int)
    seq = matrix[i,:]
    matrix[i,:] = matrix[j,:]
    matrix[j,:] = seq
    return matrix
end

function swapsequences!(matrix::NamedArray, i::Int, j::Int)
    swapsequences!(getarray(matrix), i, j)
    setnames!(matrix, _swap!(sequencenames(matrix), i, j), 1)
    return matrix
end

function swapsequences!(matrix::NamedArray, i::String, j::String)
    seqnames = sequencenames(matrix)
    swapsequences!(matrix,
                   _get_seqid_index(seqnames,i),
                   _get_seqid_index(seqnames,j))
end

"""
It puts the sequence `i` (name or position) as reference (first sequence) of the MSA. This
function swaps the sequences 1 and `i`.
"""
function setreference!(msa::AnnotatedMultipleSequenceAlignment, i::Int, annotate::Bool=true)
    swapsequences!(namedmatrix(msa), 1, i)
    if annotate
        seqnames = sequencenames(msa)
        annotate_modification!(msa, string("setreference! : Using ",
                                            seqnames[1], " instead of ",
                                            seqnames[i], " as reference."))
    end
    msa
end

function setreference!(msa::MultipleSequenceAlignment, i::Int, annotate::Bool=false)
    # The annotate argument is useful for calling this inside other functions
    swapsequences!(namedmatrix(msa), 1, i)
    msa
end

function setreference!(msa::NamedResidueMatrix{T},
                       i::Int, annotate::Bool=false) where T <: AbstractArray
    swapsequences!(msa, 1, i)
end

function setreference!(msa::NamedResidueMatrix{T},
                       id::String, annotate::Bool=false) where T <: AbstractArray
    swapsequences!(msa, names(msa,1)[1], id)
end

function setreference!(msa::AbstractMultipleSequenceAlignment, id::String,
                       annotate::Bool=true)
    setreference!(msa, _get_seqid_index(sequencenames(msa), id), annotate)
end

function setreference!(msa::Matrix{Residue}, i::Int, annotate::Bool=false)
    # The annotate argument is useful for calling this inside other functions
    swapsequences!(msa, 1, i)
end

"""
Creates a new matrix of residues. This function deletes positions/columns of the MSA with
gaps in the reference (first) sequence.
"""
function adjustreference(msa::AbstractMatrix{Residue}, annotate::Bool=false)
    # The annotate argument is useful for calling this inside other functions
    msa[:, msa[1,:] .!= GAP ]
end

"""
It removes positions/columns of the MSA with gaps in the reference (first) sequence.
"""
function adjustreference!(msa::AbstractMultipleSequenceAlignment, annotate::Bool=true)
    filtercolumns!(msa, vec(getresidues(getsequence(msa,1))) .!= GAP, annotate)
end

"""
Creates a new matrix of `Residue`s (MSA) with deleted sequences and columns/positions. The
MSA is edited in the following way:

 1. Removes all the columns/position on the MSA with gaps on the reference (first) sequence
 2. Removes all the sequences with a coverage with respect to the number of
 columns/positions on the MSA **less** than a `coveragelimit`
 (default to `0.75`: sequences with 25% of gaps)
 3. Removes all the columns/position on the MSA with **more** than a `gaplimit`
 (default to `0.5`: 50% of gaps)
"""
function gapstrip(msa::AbstractMatrix{Residue}; coveragelimit::Float64=0.75,
                  gaplimit::Float64=0.5)
    msa = adjustreference(msa)
    # Remove sequences with pour coverage of the reference sequence
    if ncolumns(msa) != 0
        msa = filtersequences(msa, vec(coverage(msa) .>= coveragelimit))
    else
        throw("There are not columns in the MSA after the gap trimming")
    end
    if nsequences(msa) != 0
        msa = filtercolumns(msa, vec(columngapfraction(msa) .<= gaplimit))
    else
        throw("There are not sequences in the MSA after coverage filter")
    end
    msa
end

"""
This functions deletes/filters sequences and columns/positions on the MSA on the following
order:

1. Removes all the columns/position on the MSA with gaps on the reference (first) sequence.
2. Removes all the sequences with a coverage with respect to the number of
columns/positions on the MSA **less** than a `coveragelimit`
(default to `0.75`: sequences with 25% of gaps).
3. Removes all the columns/position on the MSA with **more** than a `gaplimit`
(default to `0.5`: 50% of gaps).
"""
function gapstrip!(msa::AbstractMultipleSequenceAlignment,
                   annotate::Bool=isa(msa,AnnotatedAlignedObject);
                   coveragelimit::Float64=0.75, gaplimit::Float64=0.5)
    if annotate
        annotate_modification!(msa,
            string("gapstrip! : Deletes columns with gaps in the first sequence."))
    end
    adjustreference!(msa, annotate)
    # Remove sequences with pour coverage of the reference sequence
    if ncolumns(msa) != 0
        if annotate
            annotate_modification!(msa,
                string("gapstrip! : Deletes sequences with a coverage less than ",
                coveragelimit))
        end
        filtersequences!(msa, vec(coverage(msa) .>= coveragelimit), annotate)
    else
        throw("There are not columns in the MSA after the gap trimming")
    end
    # Remove columns with a porcentage of gap greater than gaplimit
    if nsequences(msa) != 0
        if annotate
            annotate_modification!(msa,
                string("gapstrip! : Deletes columns with more than ", gaplimit, " gaps."))
        end
        filtercolumns!(msa, vec(columngapfraction(msa) .<= gaplimit), annotate)
    else
        throw("There are not sequences in the MSA after coverage filter")
    end
    msa
end
