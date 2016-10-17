
# Counting Gaps and Coverage
# --------------------------

"""
Calculates the fraction of gaps on the `Array` (alignment, sequence, column, etc.).
This function can take an extra `dim` argument for calculation of the gap fraction over the given dimension
"""
function gapfraction(x::AbstractArray{Residue})
    counter = 0
    len = 0
    for res in x
        counter += ifelse(res == GAP, 1, 0)
        len += 1
    end
    float(counter) / float(len)
end

gapfraction(x::AbstractArray{Residue},
              dim::Int) = vec( mapslices(gapfraction, x, dim) )

"""
Calculates the fraction of residues (no gaps) on the `Array` (alignment, sequence, column, etc.)
This function can take an extra `dim` argument for calculation of the residue fraction over the given dimension
"""
function residuefraction(x::AbstractArray{Residue})
    counter = 0
    len = 0
    for res in x
        counter += res == GAP ? 0 : 1
        len += 1
    end
    float(counter) / float(len)
end

residuefraction(x::AbstractArray{Residue},
                  dim::Int) = vec( mapslices(residuefraction, x, dim) )

"Coverage of the sequences with respect of the number of positions on the MSA"
coverage(msa::Matrix{Residue}) = residuefraction(msa, 2)
coverage(msa::AbstractMultipleSequenceAlignment) = coverage(msa.msa)

"Fraction of gaps per column/position on the MSA"
columngapfraction(msa::Matrix{Residue}) = gapfraction(msa, 1)
columngapfraction(msa::AbstractMultipleSequenceAlignment) = columngapfraction(msa.msa)

# Reference
# ---------

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

# MSA getters
# -----------

"""
`names(msa)`

Returns the `MSA` sequences names or identifiers as an `IndexedArray`.
"""
names(msa::AbstractMultipleSequenceAlignment) = msa.id
names(msa::Matrix{Residue}) = IndexedArrays.IndexedArray(String[ string(i) for i in 1:nsequences(msa) ])

# Get annotations
# ---------------

for getter in [ :getannotcolumn, :getannotfile, :getannotresidue, :getannotsequence ]
    @eval $(getter)(msa::AnnotatedMultipleSequenceAlignment, args...) = $(getter)(msa.annotations, args...)
    @eval $(getter)(seq::AnnotatedAlignedSequence, args...) = $(getter)(seq.annotations, args...)
end

# Set annotations
# ---------------

for setter in [ :setannotcolumn!, :setannotfile!, :setannotresidue!, :setannotsequence!, :annotate_modification!, :delete_annotated_modifications!, :printmodifications ]
    @eval $(setter)(msa::AnnotatedMultipleSequenceAlignment, args...) = $(setter)(msa.annotations, args...)
    @eval $(setter)(seq::AnnotatedAlignedSequence, args...) = $(setter)(seq.annotations, args...)
end

# Used on AbstractMultipleSequenceAlignment methods
@inline annotate_modification!(msa::MultipleSequenceAlignment,
                               str::String) = false # annotate_modification! is used on bool context: annotate && ...

# Show & Print
# ------------

"""
Gives an String with the sequence number `seq` of the MSA
"""
asciisequence(msa::Matrix{Residue}, seq::Int) = ascii(convert(Vector{UInt8}, vec(msa[seq,:])))
asciisequence(msa::AbstractMultipleSequenceAlignment, seq::Int) = asciisequence(msa.msa, seq)
asciisequence(msa::AbstractMultipleSequenceAlignment, id::String) = asciisequence(msa.msa, findfirst(msa.id, id))

# Mapping annotations
# ===================

"""
Converts a string of mappings into a vector of `Int`s

```
julia> _str2int_mapping(",,2,,4,5")
6-element Array{Int64,1}:
 0
 0
 2
 0
 4
 5

```
"""
function _str2int_mapping(mapping::String)
    values = split(mapping, ',')
    len = length(values)
    intmap = Array(Int, len)
    @inbounds for i in 1:len
        value = values[i]
        intmap[i] = value == "" ? 0 : parse(Int, value)
    end
    intmap
end

getcolumnmapping(msa::AnnotatedMultipleSequenceAlignment) = _str2int_mapping(getannotfile(msa, "ColMap"))

"""
Returns the sequence coordinates as a `Vector{Int}` for an MSA sequence. That vector has one element for each MSA column.
If the number if `0` in the mapping, there is a gap in that column for that sequence.
"""
getsequencemapping(msa::AnnotatedMultipleSequenceAlignment,
                   seq_id::String) = _str2int_mapping(getannotsequence(msa, seq_id, "SeqMap"))

getsequencemapping(msa::AnnotatedMultipleSequenceAlignment,
                   seq_num::Int) = getsequencemapping(msa, msa.id[seq_num])
