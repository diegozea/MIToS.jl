# Filters
# -------

"It's similar to `filtersequences!` but for an `AbstractMatrix{Residue}`"
filtersequences(msa::AbstractMatrix{Residue}, mask::AbstractVector{Bool}) = msa[mask, :]

"""
`filtersequences!(msa, mask[, annotate::Bool=true])`

It allows to filter `msa` sequences using a `AbstractVector{Bool}` `mask`
(It removes sequnces with `false` values). `AnnotatedMultipleSequenceAlignment` annotations
are updated if `annotate` is `true` (default).
"""
function filtersequences!(msa::AnnotatedMultipleSequenceAlignment,
                          mask::AbstractVector{Bool}, annotate::Bool=true)
    msa.matrix = filtersequences(namedmatrix(msa), mask)
    filtersequences!(annotations(msa), names(msa), mask)
    annotate && annotate_modification!(msa, string("filtersequences! : ",
                                       sum(~mask), " sequences have been deleted."))
    msa
end

function filtersequences!(msa::MultipleSequenceAlignment,
                          mask::AbstractVector{Bool},
                          annotate::Bool=false) # annotate is useful for calling this
                                                # inside other functions
    msa.matrix = filtersequences(namedmatrix(msa), mask)
    msa
end

"It's similar to `filtercolumns!` but for an `AbstractMatrix{Residue}`"
filtercolumns(msa::AbstractMatrix{Residue}, mask::AbstractVector{Bool}) = msa[:, mask]

"""
`filtercolumns!(msa, mask[, annotate::Bool=true])`

It allows to filter MSA or aligned sequence columns/positions using a
`AbstractVector{Bool}` `mask`. Annotations are updated if `annotate` is `true` (default).
"""
function filtercolumns!(x::AnnotatedAlignedObject,
                        mask::AbstractVector{Bool}, annotate::Bool=true)
    x.matrix = filtercolumns(namedmatrix(x), mask)
    filtercolumns!(annotations(x), mask)
    annotate && annotate_modification!(x,string("filtercolumns! : ", sum(~mask),
                                                " columns have been deleted."))
    x
end

function filtercolumns!(x::UnannotatedAlignedObject,
                        mask::AbstractVector{Bool},
                        annotate::Bool=false)   # annotate is useful for calling this
                                                # inside other functions
    x.matrix = filtercolumns(namedmatrix(x), mask)
    x
end

filtercolumns(x::AbstractAlignedObject, args...) = filtercolumns!(deepcopy(x), args...)

# Reference
# ---------

function _swap!(names::Vector{String}, i::Int, j::Int)
    name = names[i,:]
    names[i,:] = names[j,:]
    names[j,:] = name
    names
end

"""
...
"""
function swapsequences!(matrix::NamedArray, i::Int, j::Int)
    # swap sequences
    mat = array(matrix)
    seq = mat[i,:]
    mat[i,:] = mat[j,:]
    mat[j,:] = seq
    # swap names
    setnames!(matrix, _swap!(sequencenames(matrix), i, j), 1)
    return matrix
end

function swapsequences!(matrix::NamedArray, i::String, j::String)
    seqnames = sequencenames(matrix)
    swapsequences!(matrix, findfirst(seqnames,i), findfirst(seqnames,j))
end

"""
Puts the sequence `i` as reference (as the first sequence) of the MSA.
This function swaps the sequences 1 and `i`, also an `id` can be used to select the sequence.
"""
function setreference!(msa::AnnotatedMultipleSequenceAlignment, i::Int, annotate::Bool=true)
    swap!(msa.id, 1, i)
    msa.msa[1, :], msa.msa[i, :] = msa.msa[i, :], msa.msa[1, :]
    annotate && annotate_modification!(msa, string("setreference! : Using ",
                                                   msa.id[1]," instead of ", msa.id[i]," as reference."))
    msa
end

function setreference!(msa::MultipleSequenceAlignment, i::Int, annotate::Bool=false) # annotate is useful for calling this inside other functions
    swap!(msa.id, 1, i)
    msa.msa[1, :], msa.msa[i, :] = msa.msa[i, :], msa.msa[1, :]
    msa
end

setreference!(msa::AbstractMultipleSequenceAlignment, id::String,
              annotate::Bool=true) = setreference!(msa, findfirst(msa.id ,id), annotate)

function setreference!(msa::Matrix{Residue}, i::Int, annotate::Bool=false)
    msa[1, :], msa[i, :] = msa[i, :], msa[1, :]
    msa
end

"""
Creates a new Matrix{Residue}. This function deletes positions/columns of the MSA with gaps in the reference (first) sequence
"""
adjustreference(msa::Matrix{Residue},
                annotate::Bool=false) = msa[ : , vec(msa[1,:]) .!= GAP ] # annotate is useful for calling this inside other functions

"""
Removes positions/columns of the MSA with gaps in the reference (first) sequence
"""
adjustreference!(msa::AbstractMultipleSequenceAlignment,
                 annotate::Bool=true) = filtercolumns!(msa, vec(msa.msa[1,:]) .!= GAP, annotate)

"""
This functions deletes/filters sequences and columns/positions on the MSA on the following order:

 - Removes all the columns/position on the MSA with gaps on the reference sequence (first sequence)
 - Removes all the sequences with a coverage (with respect to the number of columns/positions on the MSA) **less** than a `coveragelimit` (default to `0.75`)
 - Removes all the columns/position on the MSA with **more** than a `gaplimit` (default to `0.5`: 50% of gaps)
"""
function gapstrip!(msa::AbstractMultipleSequenceAlignment, annotate::Bool=true;
                   coveragelimit::Float64=0.75, gaplimit::Float64=0.5)
    annotate && annotate_modification!(msa, string("gapstrip! : Deletes columns with gaps in the first sequence."))
    adjustreference!(msa, annotate)
    # Remove sequences with pour coverage of the reference sequence
    if ncolumns(msa) != 0
        annotate && annotate_modification!(msa, string("gapstrip! : Deletes sequences with a coverage less than ",
                                                       coveragelimit))
        filtersequences!(msa, coverage(msa) .>= coveragelimit, annotate)
    else
        throw("There are not columns in the MSA after the gap trimming")
    end
    # Remove columns with a porcentage of gap greater than gaplimit
    if nsequences(msa) != 0
        annotate && annotate_modification!(msa, string("gapstrip! : Deletes columns with more than ",
                                                       gaplimit, " gaps."))
        filtercolumns!(msa, columngapfraction(msa) .<= gaplimit, annotate)
    else
        throw("There are not sequences in the MSA after coverage filter")
    end
    msa
end

"""
Creates a new `Matrix{Residue}` with deleted sequences and columns/positions on the MSA:

 - Removes all the columns/position on the MSA with gaps on the reference sequence (first sequence)
 - Removes all the sequences with a coverage with respect to the number of columns/positions on the MSA **less** than a `coveragelimit` (default to `0.75`: sequences with 25% of gaps)
 - Removes all the columns/position on the MSA with **more** than a `gaplimit` (default to `0.5`: 50% of gaps)
"""
function gapstrip(msa::Matrix{Residue}; coveragelimit::Float64=0.75,
                  gaplimit::Float64=0.5)
    msa = adjustreference(msa)
    # Remove sequences with pour coverage of the reference sequence
    if ncolumns(msa) != 0
        msa = filtersequences(msa, coverage(msa) .>= coveragelimit )
    else
        throw("There are not columns in the MSA after the gap trimming")
    end
    if nsequences(msa) != 0
        msa = filtercolumns(msa, columngapfraction(msa) .<= gaplimit)
    else
        throw("There are not sequences in the MSA after coverage filter")
    end
    msa
end
