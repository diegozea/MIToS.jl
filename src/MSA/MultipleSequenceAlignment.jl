import Base: length, getindex, setindex!, size, copy, deepcopy, empty!,
             convert, transpose, ctranspose, names

"""
MIToS MSAs are subtypes of `AbstractMatrix{Residue}`,
because the most basic implementation of a MIToS MSA is a `Matrix` of `Residue`s.
"""
abstract AbstractMultipleSequenceAlignment <: AbstractMatrix{Residue}

"""
MIToS sequences are subtypes of `AbstractVector{Residue}`.
"""
abstract AbstractAlignedSequence <: AbstractVector{Residue}

# Multiple Sequence Alignment
# ===========================

"""
This MSA type include the `Matrix` of `Residue`s and the sequence names.
To allow fast indexing of MSAs using **sequence identifiers**,
they are saved as an `IndexedArray`.
"""
type MultipleSequenceAlignment <: AbstractMultipleSequenceAlignment
    id::IndexedArray{ASCIIString}
    msa::Matrix{Residue}
end

"""
...
"""
type AnnotatedMultipleSequenceAlignment <: AbstractMultipleSequenceAlignment
    id::IndexedArray{ASCIIString}
    msa::Matrix{Residue}
    annotations::Annotations
end

convert(::Type{MultipleSequenceAlignment}, msa::AnnotatedMultipleSequenceAlignment) = MultipleSequenceAlignment(msa.id, msa.msa)

# Aligned Sequence
# ================

type AlignedSequence <: AbstractAlignedSequence
    id::ASCIIString
    index::Int
    sequence::Vector{Residue}
end

type AnnotatedAlignedSequence <: AbstractAlignedSequence
    id::ASCIIString
    index::Int
    sequence::Vector{Residue}
    annotations::Annotations
end

convert(::Type{AlignedSequence}, seq::AnnotatedAlignedSequence) = AlignedSequence(seq.id, seq.index, seq.sequence)

# AbstractArray Interface
# -----------------------

for meth in (:size, :length)
    @eval $(meth)(msa::AbstractMultipleSequenceAlignment) = $(meth)(msa.msa)
    @eval $(meth)(seq::AbstractAlignedSequence) = $(meth)(seq.sequence)
end

for T in (:(AlignedSequence), :(AnnotatedAlignedSequence),
          :(MultipleSequenceAlignment), :(AnnotatedMultipleSequenceAlignment))
    @eval Base.linearindexing(::Type{$(T)}) = Base.LinearFast()
end

getindex(msa::AbstractMultipleSequenceAlignment, i::Int) = getindex(msa.msa, i)
getindex(seq::AbstractAlignedSequence, i::Int) = getindex(seq.sequence, i)
setindex!(msa::AbstractMultipleSequenceAlignment, value::Residue, i::Int) =  setindex!(msa.msa, value, i)
setindex!(seq::AbstractAlignedSequence, value::Residue, i::Int) = setindex!(seq.sequence, value, i)

# Transpose
# ---------

transpose(msa::AbstractMultipleSequenceAlignment)  = transpose(msa.msa)
ctranspose(msa::AbstractMultipleSequenceAlignment) = transpose(msa.msa) # transpose is ~ 0.00022 seconds faster than ctranspose for PF00085

# Selection without Mappings
# --------------------------

"""
Allows you to access the residues in a `Matrix{Residues}`/`Vector{Residues}` without annotations.
"""
getresidues(msa::AbstractMultipleSequenceAlignment) = msa.msa
getresidues(seq::AbstractAlignedSequence) = seq.sequence

"""
Gives you the number of sequences on the `MultipleSequenceAlignment`
"""
nsequences(msa::AbstractMultipleSequenceAlignment) = size(msa.msa, 1)
nsequences(msa::Matrix{Residue}) = size(msa, 1)

"""
Gives you the number of columns/positions on the MSA or aligned sequence
"""
ncolumns(msa::AbstractMultipleSequenceAlignment) = size(msa.msa, 2)
ncolumns(msa::Matrix{Residue}) = size(msa, 2)
ncolumns(seq::AbstractAlignedSequence) = length(seq.sequence)
ncolumns(seq::Vector{Residue}) = length(seq)

"""
Gives you a `Vector{Vector{Residue}}` with all the sequences of the MSA without Annotations
"""
function getresiduesequences(msa::Matrix{Residue})
    nseq = nsequences(msa)
    tmsa = msa'
    sequences = Array(Vector{Residue}, nseq)
    for i in 1:nseq
        @inbounds sequences[i] = tmsa[:,i]
    end
    sequences
end

getresiduesequences(msa::AbstractMultipleSequenceAlignment) = getresiduesequences(msa.msa)

# Select sequence
# ---------------
"""
Gives you the annotations of the Sequence
"""
function getsequence(data::Annotations, id::ASCIIString)
    GS = Dict{Tuple{ASCIIString,ASCIIString},ASCIIString}()
    GR = Dict{Tuple{ASCIIString,ASCIIString},ASCIIString}()
    if length(data.sequences) > 0 || length(data.residues) > 0
        for (key, value) in data.sequences
            if key[1] == id
                GS[key] = value
            end
        end
        for (key, value) in data.residues
            if key[1] == id
                GR[key] = value
            end
        end
        GS = sizehint!(GS, length(GS))
        GR = sizehint!(GR, length(GR))
    end
    Annotations(data.file, GS, data.columns, GR)
end

"""
Returns an `AnnotatedAlignedSequence` with all annotations of sequence from the `AnnotatedMultipleSequenceAlignment`
"""
getsequence(msa::AnnotatedMultipleSequenceAlignment,
            i::Int) = AnnotatedAlignedSequence(msa.id[i], i,  vec(msa.msa[i,:]),
                                               getsequence(msa.annotations, msa.id[i]))

"""
Returns an `AlignedSequence` from the `MultipleSequenceAlignment`
"""
getsequence(msa::MultipleSequenceAlignment,
            i::Int) = AlignedSequence(msa.id[i], i, vec(msa.msa[i,:]))

getsequence(msa::Matrix{Residue}, i::Int) = vec(msa[i,:])

function getsequence(msa::AnnotatedMultipleSequenceAlignment, id::ASCIIString)
    i = findfirst(msa.id, id)
    AnnotatedAlignedSequence(id, i, vec(msa.msa[i,:]),
                             getsequence(msa.annotations, id))
end

function getsequence(msa::MultipleSequenceAlignment, id::ASCIIString)
    i = findfirst(msa.id, id)
    AlignedSequence(id, i, vec(msa.msa[i,:]))
end

getindex(msa::AbstractMultipleSequenceAlignment, id::ASCIIString) = getsequence(msa, id)

# Filters
# -------

filtersequences(msa::Matrix{Residue}, mask::AbstractVector{Bool}) = msa[mask, :]

"""
Allows to filter sequences on a MSA using a `AbstractVector{Bool}` mask (removes `false`s).
For `AnnotatedMultipleSequenceAlignment`s the annotations are updated.
"""
function filtersequences!(msa::AnnotatedMultipleSequenceAlignment,
                          mask::AbstractVector{Bool}, annotate::Bool=true)
    msa.msa = filtersequences(msa.msa, mask)
    #msa.sequencemapping = msa.sequencemapping[ mask , : ]
    filtersequences!(msa.annotations, msa.id, mask)
    msa.id = IndexedArray(msa.id[ mask ])
    annotate && annotate_modification!(msa, string("filtersequences! : ", sum(~mask), " sequences have been deleted."))
    msa
end

function filtersequences!(msa::MultipleSequenceAlignment,
                          mask::AbstractVector{Bool}, annotate::Bool=false) # annotate is useful for calling this inside other functions
    msa.msa = filtersequences(msa.msa, mask)
    msa.id = IndexedArray(msa.id[ mask ])
    msa
end

filtercolumns(msa::Matrix{Residue}, mask::AbstractVector{Bool}) = msa[ : , mask ]
filtercolumns(seq::Vector{Residue}, mask::AbstractVector{Bool}) = seq[ mask ]

"""
Allows to filter columns/positions on a MSA using a `AbstractVector{Bool}` mask.
For `AnnotatedMultipleSequenceAlignment`s or `AnnotatedAlignedSequence`s the annotations are updated.
"""
function filtercolumns!(msa::AnnotatedMultipleSequenceAlignment,
                        mask::AbstractVector{Bool}, annotate::Bool=true)
    msa.msa = filtercolumns(msa.msa, mask)
    #msa.sequencemapping = msa.sequencemapping[ : , mask ]
    #msa.filecolumnmapping = msa.filecolumnmapping[ mask ]
    filtercolumns!(msa.annotations, mask)
    annotate && annotate_modification!(msa, string("filtercolumns! : ",
                                                   sum(~mask), " columns have been deleted."))
    msa
end

function filtercolumns!(msa::MultipleSequenceAlignment, mask::AbstractVector{Bool}, annotate::Bool=false) # annotate is useful for calling this inside other functions
    msa.msa = filtercolumns(msa.msa, mask)
    msa
end

function filtercolumns!(seq::AnnotatedAlignedSequence, mask::AbstractVector{Bool}, annotate::Bool=true)
    seq.sequence = filtercolumns(seq.sequence, mask)
    #seq.sequencemapping = seq.sequencemapping[ mask ]
    #seq.filecolumnmapping = seq.filecolumnmapping[ mask ]
    filtercolumns!(seq.annotations, mask)
    annotate && annotate_modification!(seq, string("filtercolumns! : ",
                                                   sum(~mask), " columns have been deleted."))
    seq
end

function filtercolumns!(seq::AlignedSequence, mask::AbstractVector{Bool}, annotate::Bool=false)
    seq.sequence = filtercolumns(seq.sequence, mask)
    seq
end

# Copy, deepcopy, empty!
# ----------------------

for fun in (:copy, :deepcopy)
    @eval $(fun)(msa::AnnotatedMultipleSequenceAlignment) = AnnotatedMultipleSequenceAlignment($(fun)(msa.id), $(fun)(msa.msa), $(fun)(msa.annotations))
    @eval $(fun)(msa::MultipleSequenceAlignment) = MultipleSequenceAlignment($(fun)(msa.id), $(fun)(msa.msa))
    @eval $(fun)(seq::AnnotatedAlignedSequence) = AnnotatedAlignedSequence($(fun)(seq.id), $(fun)(seq.index), $(fun)(seq.sequence), $(fun)(seq.annotations))
    @eval $(fun)(seq::AlignedSequence) = AlignedSequence($(fun)(seq.id), $(fun)(seq.index), $(fun)(seq.sequence))
end

empty!(msa::AnnotatedMultipleSequenceAlignment) = (empty!(msa.id); empty!(msa.msa);
                                                   empty!(msa.annotations); msa)
empty!(msa::MultipleSequenceAlignment) = (empty!(msa.id); empty!(msa.msa); msa)
empty!(seq::AnnotatedAlignedSequence) = (empty!(seq.id); empty!(seq.index);
                                         v(seq.sequence); empty!(seq.annotations); seq)
empty!(seq::AlignedSequence) = (empty!(seq.id); empty!(seq.index);
                                empty!(seq.sequence); seq)

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
        counter += res == GAP ? 1 : 0
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

setreference!(msa::AbstractMultipleSequenceAlignment, id::ASCIIString,
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

"Returns the annotations of a MSA or a sequence."
annotations(msa::AnnotatedMultipleSequenceAlignment) = msa.annotations
annotations(seq::AnnotatedAlignedSequence) = seq.annotations

"Returns the names of the MSA sequences."
names(msa::AbstractMultipleSequenceAlignment) = msa.id

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
                               str::ASCIIString) = false # annotate_modification! is used on bool context: annotate && ...

# Show & Print
# ------------

"""
Gives an ASCIIString with the sequence number `seq` of the MSA
"""
asciisequence(msa::Matrix{Residue}, seq::Int) = ascii(convert(Vector{UInt8}, vec(msa[seq,:])))
asciisequence(msa::AbstractMultipleSequenceAlignment, seq::Int) = asciisequence(msa.msa, seq)

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
function _str2int_mapping(mapping::ASCIIString)
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

getsequencemapping(msa::AnnotatedMultipleSequenceAlignment,
                   seq_id::ASCIIString) = _str2int_mapping(getannotsequence(msa, seq_id, "SeqMap"))

getsequencemapping(msa::AnnotatedMultipleSequenceAlignment,
                   seq_num::Int) = getsequencemapping(msa, msa.id[seq_num])
