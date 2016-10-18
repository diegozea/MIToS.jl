
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
