import Base: length, getindex, size, ndims

# Multiple Sequence Alignment
# ===========================

type MultipleSequenceAlignment
  id::IndexedVector{ASCIIString}
  msa::Matrix{Residue}
  sequencemapping::Matrix{Int}
  filecolumnmapping::IndexedVector{Int}
  annotations::Annotations
end

# Aligned Sequence
# ================

type AlignedSequence
  id::ASCIIString
  index::Int
  sequence::Vector{Residue}
  sequencemapping::Vector{Int}
  filecolumnmapping::IndexedVector{Int}
  annotations::Annotations
end

# Selection without Mappings
# --------------------------

getresidues(msa::MultipleSequenceAlignment) = msa.msa

getresidues(msa::MultipleSequenceAlignment, sequence::Int, column::Int) = msa.msa[sequence, column]

getresidues(msa::MultipleSequenceAlignment, sequence::Int, range::UnitRange) = msa.msa[sequence, range]
getresidues(msa::MultipleSequenceAlignment, sequence::Int, range::Colon) = msa.msa[sequence, :]

getresidues(msa::MultipleSequenceAlignment, range::UnitRange, column::Int) = msa.msa[range, column]
getresidues(msa::MultipleSequenceAlignment, range::Colon, column::Int) = msa.msa[:, column]

size(msa::MultipleSequenceAlignment, i::Int) = size(msa.msa, i)
size(msa::MultipleSequenceAlignment) = size(msa.msa)

ndims(msa::MultipleSequenceAlignment) = ndims(msa.msa)

nsequences(msa::MultipleSequenceAlignment) = size(msa, 1)
ncolumns(msa::MultipleSequenceAlignment) = size(msa, 2)

# Select sequence
# ---------------

function getsequence(data::Annotations, id::ASCIIString)
  GS = Dict{(ASCIIString,ASCIIString),ASCIIString}()
  GR = Dict{(ASCIIString,ASCIIString),ASCIIString}()
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
    GS = sizehint(GS, length(GS))
    GR = sizehint(GR, length(GR))
  end
  Annotations(data.file, GS, data.columns, GR)
end

function getsequence(msa::MultipleSequenceAlignment, id::ASCIIString)
  i = selectindex(msa.id, id)
  AlignedSequence(id, i, vec(msa.msa[i,:]), vec(msa.sequencemapping[i,:]), msa.filecolumnmapping, getsequence(msa.annotations, id))
end

getsequence(msa::MultipleSequenceAlignment, i::Int) = AlignedSequence(selectvalue(msa.id, i), i,  vec(msa.msa[i,:]), vec(msa.sequencemapping[i,:]), msa.filecolumnmapping, getsequence(msa.annotations, selectvalue(msa.id, i)))

getindex(msa::MultipleSequenceAlignment, id::ASCIIString) = getsequence(msa, id)

# Filters
# -------

function filtersequences!(msa::MultipleSequenceAlignment, mask::BitVector)
  if length(mask) == nsequences(msa)
    msa.id = msa.id[ mask ]
    msa.msa = msa.msa[ mask , : ]
    msa.sequencemapping = msa.sequencemapping[ mask , : ]
    filtersequences!(msa.annotations, msa.id, mask)
  end
  msa
end

function filtercolumns!(msa::MultipleSequenceAlignment, mask::BitVector)
  if length(mask) == ncolumns(msa)
    msa.msa = msa.msa[ : , mask ]
    msa.sequencemapping = msa.sequencemapping[ : , mask ]
    msa.filecolumnmapping = msa.filecolumnmapping[ mask ]
    filtercolumns!(msa.annotations, mask)
  end
  msa
end

# Copy and deepcopy
# -----------------

deepcopy(msa::MultipleSequenceAlignment) = MultipleSequenceAlignment( deepcopy( msa.id ), deepcopy( msa.msa ),
                                                                     deepcopy( msa.sequencemapping ), deepcopy( msa.filecolumnmapping ), deepcopy(msa.annotations) )

copy(msa::MultipleSequenceAlignment) = MultipleSequenceAlignment( copy( msa.id ), copy( msa.msa ),
                                                                 copy( msa.sequencemapping ), copy( msa.filecolumnmapping ), copy(msa.annotations))

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

"""Coverage of sequences"""
coverage(msa::MultipleSequenceAlignment) = vec( mapslices(gappercentage, msa.msa, 2) )

columngappercentage(msa::MultipleSequenceAlignment) = vec( mapslices(gappercentage, msa.msa, 1) )

# Reference
# ---------

function setreference!(msa::MultipleSequenceAlignment, i::Int)
  swap!(msa.id, 1, i)
  msa.msa[1, :], msa.msa[i, :] = msa.msa[i, :], msa.msa[1, :]
  msa.sequencemapping[1, :], msa.sequencemapping[i, :] = msa.sequencemapping[i, :], msa.sequencemapping[1, :]
  msa
end

setreference!(msa::MultipleSequenceAlignment, id::ASCIIString) = setreference!(msa, selectindex(msa.id ,id))

function gapstrip!(msa::MultipleSequenceAlignment; coveragelimit::FloatingPoint=0.8, gaplimit::FloatingPoint=0.5)
  adjustreference!(msa)
  # Remove sequences with pour coverage of the reference sequence
  ncolumns(msa) != 0 ? filtersequences!(msa, coverage(msa) .>= coveragelimit ) : throw("There are not columns in the MSA after the gap trimming")
  # Remove columns with a porcentage of gap greater than gaplimit
  nsequences(msa) != 0 ? filtercolumns!(msa, columngappercentage(msa) .<= gaplimit) : throw("There are not sequences in the MSA after coverage filter")
end

# Remove positions with gaps in the reference sequence
adjustreference!(msa::MultipleSequenceAlignment) = filtercolumns!(msa, vec(msa.msa[1,:]) .!= GAP )
