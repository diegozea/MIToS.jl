import Base: length, getindex, setindex!, size, copy, deepcopy, empty!, convert

abstract AbstractMultipleSequenceAlignment <: AbstractArray{Residue,2}
abstract AbstractSequence <: AbstractArray{Residue,1}

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

type Sequence <: AbstractSequence
  id::ASCIIString
  index::Int
  sequence::Vector{Residue}
end

type AnnotatedSequence <: AbstractSequence
  id::ASCIIString
  index::Int
  sequence::Vector{Residue}
  sequencemapping::Vector{Int}
  filecolumnmapping::IndexedVector{Int}
  annotations::Annotations
end

convert(::Type{Sequence}, seq::AnnotatedSequence) = Sequence(seq.id, seq.index, seq.sequence)

# AbstractArray Interface
# -----------------------

for meth in (:size, :length)
  @eval $(meth)(msa::AbstractMultipleSequenceAlignment) = $(meth)(msa.msa)
  @eval $(meth)(seq::AbstractSequence) = $(meth)(seq.sequence)
end

for T in (:(Sequence), :(AnnotatedSequence), :(MultipleSequenceAlignment), :(AnnotatedMultipleSequenceAlignment))
  @eval Base.linearindexing(::Type{$(T)}) = Base.LinearFast()
end

getindex(msa::AbstractMultipleSequenceAlignment, i::Int) = getindex(msa.msa, i)
getindex(seq::AbstractSequence, i::Int) = getindex(seq.sequence, i)
setindex!(msa::AbstractMultipleSequenceAlignment, value::Residue, i::Int) =  setindex!(msa.msa, value, i)
setindex!(seq::AbstractSequence, value::Residue, i::Int) = setindex!(seq.sequence, value, i)

#getresidues(msa::MultipleSequenceAlignment, linear::Int) = getindex(msa.msa, linear)
#getresidues(msa::MultipleSequenceAlignment, sequence::Int, column::Int) = getindex(msa.msa, sequence, column)
#getresidues(msa::MultipleSequenceAlignment, indexes...) = getindex(msa.msa, indexes...)

# size(msa::MultipleSequenceAlignment, i::Int) = size(msa.msa, i)
# ndims(msa::MultipleSequenceAlignment) = ndims(msa.msa)

# Selection without Mappings
# --------------------------

"""Allows you to access the residues in a `Matrix{Residues}`/`Vector{Residues}` without annotations."""
getresidues(msa::AbstractMultipleSequenceAlignment) = msa.msa
getresidues(seq::AbstractSequence) = seq.sequence

"""Gives you the number of sequences on the `MultipleSequenceAlignment`"""
nsequences(msa::AbstractMultipleSequenceAlignment) = size(msa, 1)

"""Gives you the number of columns/positions on the MSA or aligned sequence"""
ncolumns(msa::AbstractMultipleSequenceAlignment) = size(msa, 2)
ncolumns(seq::AbstractSequence) = length(seq.sequence)

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

"""Returns a `Sequence` or `AnnotatedSequence` with all annotations of sequence from the MSA"""
getsequence(msa::AnnotatedMultipleSequenceAlignment, i::Int) = AnnotatedSequence(selectvalue(msa.id, i), i,  vec(msa.msa[i,:]), vec(msa.sequencemapping[i,:]), msa.filecolumnmapping, getsequence(msa.annotations, selectvalue(msa.id, i)))
getsequence(msa::MultipleSequenceAlignment, i::Int) = Sequence(selectvalue(msa.id, i), i,  vec(msa.msa[i,:]))

function getsequence(msa::AnnotatedMultipleSequenceAlignment, id::ASCIIString)
  i = selectindex(msa.id, id)
  AnnotatedSequence(id, i, vec(msa.msa[i,:]), vec(msa.sequencemapping[i,:]), msa.filecolumnmapping, getsequence(msa.annotations, id))
end

function getsequence(msa::MultipleSequenceAlignment, id::ASCIIString)
  i = selectindex(msa.id, id)
  Sequence(id, i, vec(msa.msa[i,:]))
end

getindex(msa::AbstractMultipleSequenceAlignment, id::ASCIIString) = getsequence(msa, id)

# Filters
# -------

function filtersequences!(msa::AnnotatedMultipleSequenceAlignment, mask::BitVector)
  if length(mask) == nsequences(msa)
    msa.msa = msa.msa[ mask , : ]
    msa.sequencemapping = msa.sequencemapping[ mask , : ]
    filtersequences!(msa.annotations, msa.id, mask)
    msa.id = msa.id[ mask ]
    return(msa)
  else
    throw(ArgumentError("The number of sequences on the MSA should be equal to the length of the mask"))
  end
end

function filtersequences!(msa::MultipleSequenceAlignment, mask::BitVector)
  if length(mask) == nsequences(msa)
    msa.msa = msa.msa[ mask , : ]
    msa.id = msa.id[ mask ]
    return(msa)
  else
    throw(ArgumentError("The number of sequences on the MSA should be equal to the length of the mask"))
  end
end

function filtercolumns!(msa::AnnotatedMultipleSequenceAlignment, mask::BitVector)
  if length(mask) == ncolumns(msa)
    msa.msa = msa.msa[ : , mask ]
    msa.sequencemapping = msa.sequencemapping[ : , mask ]
    msa.filecolumnmapping = msa.filecolumnmapping[ mask ]
    filtercolumns!(msa.annotations, mask)
    return(msa)
  else
    throw(ArgumentError("The number of columns on the MSA should be equal to the length of the mask"))
  end
end

function filtercolumns!(msa::MultipleSequenceAlignment, mask::BitVector)
  if length(mask) == ncolumns(msa)
    msa.msa = msa.msa[ : , mask ]
    return(msa)
  else
    throw(ArgumentError("The number of columns on the MSA should be equal to the length of the mask"))
  end
end

function filtercolumns!(seq::AnnotatedSequence, mask::BitVector)
  if length(mask) == length(seq)
    seq.sequence = seq.sequence[ : , mask ]
    seq.sequencemapping = seq.sequencemapping[ : , mask ]
    seq.filecolumnmapping = seq.filecolumnmapping[ mask ]
    filtercolumns!(seq.annotations, mask)
    return(seq)
  else
    throw(ArgumentError("The number of columns/residues on the aligned sequence should be equal to the length of the mask"))
  end
end

function filtercolumns!(seq::Sequence, mask::BitVector)
  if length(mask) == length(seq)
    seq.sequence = seq.sequence[ : , mask ]
    return(seq)
  else
    throw(ArgumentError("The number of columns/residues on the aligned sequence should be equal to the length of the mask"))
  end
end

# Copy, deepcopy, empty!
# ----------------------

for fun in (:copy, :deepcopy, :empty!)
  @eval $(fun)(msa::AnnotatedMultipleSequenceAlignment) = AnnotatedMultipleSequenceAlignment($(fun)(msa.id), $(fun)(msa.msa), $(fun)(msa.sequencemapping), $(fun)(msa.filecolumnmapping), $(fun)(msa.annotations))
  @eval $(fun)(msa::MultipleSequenceAlignment) = MultipleSequenceAlignment($(fun)(msa.id), $(fun)(msa.msa))
  @eval $(fun)(seq::AnnotatedSequence) = AnnotatedSequence($(fun)(seq.id), $(fun)(seq.index), $(fun)(seq.sequence), $(fun)(seq.sequencemapping), $(fun)(seq.filecolumnmapping), $(fun)(seq.annotations))
  @eval $(fun)(seq::Sequence) = Sequence($(fun)(seq.id), $(fun)(seq.index), $(fun)(seq.sequence))
end

# Counting Gaps and Coverage
# --------------------------

function gappercentage(x::AbstractArray{Residue})
  counter = 0
  len = 0
  for res in x
    counter += res == GAP ? 1 : 0
    len += 1
  end
  float(counter) / float(len)
end

function residuepercentage(x::AbstractArray{Residue})
  counter = 0
  len = 0
  for res in x
    counter += res == GAP ? 0 : 1
    len += 1
  end
  float(counter) / float(len)
end

"""Coverage of sequences"""
coverage(msa::AbstractMultipleSequenceAlignment) = vec( mapslices(residuepercentage, msa.msa, 2) )

columngappercentage(msa::AbstractMultipleSequenceAlignment) = vec( mapslices(gappercentage, msa.msa, 1) )

# Reference
# ---------

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

function gapstrip!(msa::AbstractMultipleSequenceAlignment; coveragelimit::Float64=0.75, gaplimit::Float64=0.5)
  adjustreference!(msa)
  # Remove sequences with pour coverage of the reference sequence
  ncolumns(msa) != 0 ? filtersequences!(msa, coverage(msa) .>= coveragelimit ) : throw("There are not columns in the MSA after the gap trimming")
  # Remove columns with a porcentage of gap greater than gaplimit
  nsequences(msa) != 0 ? filtercolumns!(msa, columngappercentage(msa) .<= gaplimit) : throw("There are not sequences in the MSA after coverage filter")
end

# Remove positions with gaps in the reference sequence
adjustreference!(msa::AbstractMultipleSequenceAlignment) = filtercolumns!(msa, vec(msa.msa[1,:]) .!= GAP )

# Show & Print
# ------------

asciisequence(msa::AbstractMultipleSequenceAlignment, seq::Int) = ascii(convert(Vector{UInt8}, vec(msa.msa[seq,:])))

# print(io::IO, msa::AbstractMultipleSequenceAlignment) = dump(io, msa)
# show(io::IO, msa::AbstractMultipleSequenceAlignment) = dump(io, msa)
