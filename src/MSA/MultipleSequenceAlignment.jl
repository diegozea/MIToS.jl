import Base: length, getindex, setindex!, size, copy, deepcopy, empty!, convert

abstract AbstractMultipleSequenceAlignment <: AbstractArray{Residue,2}
abstract AbstractAlignedSequence <: AbstractArray{Residue,1}

# Multiple Sequence Alignment
# ===========================

type MultipleSequenceAlignment <: AbstractMultipleSequenceAlignment
  id::IndexedVector{ASCIIString}
  msa::Matrix{Residue}
end

type AnnotatedMultipleSequenceAlignment <: AbstractMultipleSequenceAlignment
  id::IndexedVector{ASCIIString}
  msa::Matrix{Residue}
  sequencemapping::Matrix{Int}
  filecolumnmapping::IndexedVector{Int}
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
  sequencemapping::Vector{Int}
  filecolumnmapping::IndexedVector{Int}
  annotations::Annotations
end

convert(::Type{AlignedSequence}, seq::AnnotatedAlignedSequence) = AlignedSequence(seq.id, seq.index, seq.sequence)

# AbstractArray Interface
# -----------------------

for meth in (:size, :length)
  @eval $(meth)(msa::AbstractMultipleSequenceAlignment) = $(meth)(msa.msa)
  @eval $(meth)(seq::AbstractAlignedSequence) = $(meth)(seq.sequence)
end

for T in (:(AlignedSequence), :(AnnotatedAlignedSequence), :(MultipleSequenceAlignment), :(AnnotatedMultipleSequenceAlignment))
  @eval Base.linearindexing(::Type{$(T)}) = Base.LinearFast()
end

getindex(msa::AbstractMultipleSequenceAlignment, i::Int) = getindex(msa.msa, i)
getindex(seq::AbstractAlignedSequence, i::Int) = getindex(seq.sequence, i)
setindex!(msa::AbstractMultipleSequenceAlignment, value::Residue, i::Int) =  setindex!(msa.msa, value, i)
setindex!(seq::AbstractAlignedSequence, value::Residue, i::Int) = setindex!(seq.sequence, value, i)

# Selection without Mappings
# --------------------------

"""Allows you to access the residues in a `Matrix{Residues}`/`Vector{Residues}` without annotations."""
getresidues(msa::AbstractMultipleSequenceAlignment) = msa.msa
getresidues(seq::AbstractAlignedSequence) = seq.sequence

"""Gives you the number of sequences on the `MultipleSequenceAlignment`"""
nsequences(msa::AbstractMultipleSequenceAlignment) = size(msa.msa, 1)
nsequences(msa::Matrix{Residue}) = size(msa, 1)

"""Gives you the number of columns/positions on the MSA or aligned sequence"""
ncolumns(msa::AbstractMultipleSequenceAlignment) = size(msa.msa, 2)
ncolumns(msa::Matrix{Residue}) = size(msa, 2)
ncolumns(seq::AbstractAlignedSequence) = length(seq.sequence)
ncolumns(seq::Vector{Residue}) = length(seq)

"""Gives you a `Vector{Vector{Residue}}` with all the sequences of the MSA without Annotations"""
function getresiduesequences(msa::AbstractMultipleSequenceAlignment)
  nseq = nsequences(msa)
  tmsa = msa.msa'
  sequences = Array(Vector{Residue}, nseq)
  for i in 1:nseq
    @inbounds sequences[i] = tmsa[:,i]
  end
  sequences
end

# Select sequence
# ---------------
"""Gives you the annotations of the Sequence"""
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

"""Returns an `AnnotatedAlignedSequence` with all annotations of sequence from the `AnnotatedMultipleSequenceAlignment`"""
getsequence(msa::AnnotatedMultipleSequenceAlignment, i::Int) = AnnotatedAlignedSequence(selectvalue(msa.id, i), i,  vec(msa.msa[i,:]),
                                                                                        vec(msa.sequencemapping[i,:]), msa.filecolumnmapping,
                                                                                        getsequence(msa.annotations, selectvalue(msa.id, i)))

"""Returns an `AlignedSequence` from the `MultipleSequenceAlignment`"""
getsequence(msa::MultipleSequenceAlignment, i::Int) = AlignedSequence(selectvalue(msa.id, i), i,  vec(msa.msa[i,:]))

getsequence(msa::Matrix{Residue}, i::Int) = vec(msa[i,:])

function getsequence(msa::AnnotatedMultipleSequenceAlignment, id::ASCIIString)
  i = selectindex(msa.id, id)
  AnnotatedAlignedSequence(id, i, vec(msa.msa[i,:]), vec(msa.sequencemapping[i,:]), msa.filecolumnmapping, getsequence(msa.annotations, id))
end

function getsequence(msa::MultipleSequenceAlignment, id::ASCIIString)
  i = selectindex(msa.id, id)
  AlignedSequence(id, i, vec(msa.msa[i,:]))
end

getindex(msa::AbstractMultipleSequenceAlignment, id::ASCIIString) = getsequence(msa, id)

# Filters
# -------

function filtersequences(msa::Matrix{Residue}, mask::BitVector)
  if length(mask) == nsequences(msa)
    return(msa[mask, :])
  else
    throw(ArgumentError("The number of sequences on the MSA should be equal to the length of the mask"))
  end
end

"""
Allows to filter sequences on a MSA using a `BitVector` mask (removes `false`s).
For `AnnotatedMultipleSequenceAlignment`s the annotations are updated.
"""
function filtersequences!(msa::AnnotatedMultipleSequenceAlignment, mask::BitVector)
  msa.msa = filtersequences(msa.msa, mask)
  msa.sequencemapping = msa.sequencemapping[ mask , : ]
  filtersequences!(msa.annotations, msa.id, mask)
  msa.id = msa.id[ mask ]
  msa
end

function filtersequences!(msa::MultipleSequenceAlignment, mask::BitVector)
  msa.msa = filtersequences(msa.msa, mask)
  msa.id = msa.id[ mask ]
  msa
end

function filtercolumns(msa::Matrix{Residue}, mask::BitVector)
  if length(mask) == ncolumns(msa)
    return(msa[ : , mask ])
  else
    throw(ArgumentError("The number of columns on the MSA should be equal to the length of the mask"))
  end
end

function filtercolumns(seq::Vector{Residue}, mask::BitVector)
  if length(mask) == ncolumns(seq)
    return(seq[ mask ])
  else
    throw(ArgumentError("The number of columns on the MSA should be equal to the length of the mask"))
  end
end

"""
Allows to filter columns/positions on a MSA using a `BitVector` mask.
For `AnnotatedMultipleSequenceAlignment`s or `AnnotatedAlignedSequence`s the annotations are updated.
"""
function filtercolumns!(msa::AnnotatedMultipleSequenceAlignment, mask::BitVector)
  msa.msa = filtercolumns(msa.msa, mask)
  msa.sequencemapping = msa.sequencemapping[ : , mask ]
  msa.filecolumnmapping = msa.filecolumnmapping[ mask ]
  filtercolumns!(msa.annotations, mask)
  msa
end

function filtercolumns!(msa::MultipleSequenceAlignment, mask::BitVector)
  msa.msa = filtercolumns(msa.msa, mask)
  msa
end

function filtercolumns!(seq::AnnotatedAlignedSequence, mask::BitVector)
  seq.sequence = filtercolumns(seq.sequence, mask)
  seq.sequencemapping = seq.sequencemapping[ mask ]
  seq.filecolumnmapping = seq.filecolumnmapping[ mask ]
  filtercolumns!(seq.annotations, mask)
  seq
end

function filtercolumns!(seq::AlignedSequence, mask::BitVector)
  seq.sequence = filtercolumns(seq.sequence, mask)
  seq
end

# Copy, deepcopy, empty!
# ----------------------

for fun in (:copy, :deepcopy, :empty!)
  @eval $(fun)(msa::AnnotatedMultipleSequenceAlignment) = AnnotatedMultipleSequenceAlignment($(fun)(msa.id), $(fun)(msa.msa), $(fun)(msa.sequencemapping), $(fun)(msa.filecolumnmapping), $(fun)(msa.annotations))
  @eval $(fun)(msa::MultipleSequenceAlignment) = MultipleSequenceAlignment($(fun)(msa.id), $(fun)(msa.msa))
  @eval $(fun)(seq::AnnotatedAlignedSequence) = AnnotatedAlignedSequence($(fun)(seq.id), $(fun)(seq.index), $(fun)(seq.sequence), $(fun)(seq.sequencemapping), $(fun)(seq.filecolumnmapping), $(fun)(seq.annotations))
  @eval $(fun)(seq::AlignedSequence) = AlignedSequence($(fun)(seq.id), $(fun)(seq.index), $(fun)(seq.sequence))
end

# Counting Gaps and Coverage
# --------------------------

"""Calculates the percentage of gaps on the `Array` (alignment, sequence, column, etc.)"""
function gappercentage(x::AbstractArray{Residue})
  counter = 0
  len = 0
  for res in x
    counter += res == GAP ? 1 : 0
    len += 1
  end
  float(counter) / float(len)
end

"""Calculates the percentage of residues (no gaps) on the `Array` (alignment, sequence, column, etc.)"""
function residuepercentage(x::AbstractArray{Residue})
  counter = 0
  len = 0
  for res in x
    counter += res == GAP ? 0 : 1
    len += 1
  end
  float(counter) / float(len)
end

"""Coverage of the sequences with respect of the number of positions on the MSA"""
coverage(msa::Matrix{Residue}) = vec( mapslices(residuepercentage, msa, 2) )
coverage(msa::AbstractMultipleSequenceAlignment) = coverage(msa.msa)

"""Percentage of gaps per column/position on the MSA"""
columngappercentage(msa::Matrix{Residue}) = vec( mapslices(gappercentage, msa, 1) )
columngappercentage(msa::AbstractMultipleSequenceAlignment) = columngappercentage(msa.msa)

# Reference
# ---------

"""
Puts the sequence `i` as reference (as the first sequence) of the MSA.
This function swaps the sequences 1 and `i`, also and `id` can be used.
"""
function setreference!(msa::AnnotatedMultipleSequenceAlignment, i::Int)
  swap!(msa.id, 1, i)
  msa.msa[1, :], msa.msa[i, :] = msa.msa[i, :], msa.msa[1, :]
  msa.sequencemapping[1, :], msa.sequencemapping[i, :] = msa.sequencemapping[i, :], msa.sequencemapping[1, :]
  msa
end

function setreference!(msa::MultipleSequenceAlignment, i::Int)
  swap!(msa.id, 1, i)
  msa.msa[1, :], msa.msa[i, :] = msa.msa[i, :], msa.msa[1, :]
  msa
end

setreference!(msa::AbstractMultipleSequenceAlignment, id::ASCIIString) = setreference!(msa, selectindex(msa.id ,id))

function setreference!(msa::Matrix{Residue}, i::Int)
  msa[1, :], msa[i, :] = msa[i, :], msa[1, :]
  msa
end

adjustreference(msa::Matrix{Residue}) = msa[ vec(msa[1,:]) .!= GAP ]

"""Removes positions/columns of the MSA with gaps in the reference (first) sequence"""
adjustreference!(msa::AbstractMultipleSequenceAlignment) = filtercolumns!(msa, vec(msa.msa[1,:]) .!= GAP )


"""
This functions deletes/filters sequences and columns/positions on the MSA on the following order:

 - Removes all the columns/position on the MSA with gaps on the reference sequence (first sequence)
 - Removes all the sequences with a coverage with respect to the number of columns/positions on the MSA **less** than a `coveragelimit` (default to `0.75`: sequences with 25% of gaps)
 - Removes all the columns/position on the MSA with **more** than a `gaplimit` (default to `0.5`: 50% of gaps)
"""
function gapstrip!(msa::AbstractMultipleSequenceAlignment; coveragelimit::Float64=0.75, gaplimit::Float64=0.5)
  adjustreference!(msa)
  # Remove sequences with pour coverage of the reference sequence
  ncolumns(msa) != 0 ? filtersequences!(msa, coverage(msa) .>= coveragelimit ) : throw("There are not columns in the MSA after the gap trimming")
  # Remove columns with a porcentage of gap greater than gaplimit
  nsequences(msa) != 0 ? filtercolumns!(msa, columngappercentage(msa) .<= gaplimit) : throw("There are not sequences in the MSA after coverage filter")
end

# Get annotations
# ---------------

for getter in [ :getannotcolumn, :getannotfile, :getannotresidue, :getannotsequence ]
  @eval $(getter)(msa::AnnotatedMultipleSequenceAlignment, args...) = $(getter)(msa.annotations, args...)
  @eval $(getter)(seq::AnnotatedAlignedSequence, args...) = $(getter)(seq.annotations, args...)
end

# Show & Print
# ------------

asciisequence(msa::AbstractMultipleSequenceAlignment, seq::Int) = ascii(convert(Vector{UInt8}, vec(msa.msa[seq,:])))
