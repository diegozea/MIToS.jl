# import Base: length, getindex, setindex!, size, copy, deepcopy, empty!,
#              convert, transpose, ctranspose, names

"""
MIToS MSA and aligned sequences (aligned objects) are subtypes of `AbstractMatrix{Residue}`,
because MSAs and sequences are stored as `Matrix` of `Residue`s.
"""
abstract type AbstractAlignedObject <: AbstractMatrix{Residue} end

"""
MSAs are stored as `Matrix{Residue}`. It's possible to use a `NamedArray{Residue,2}` as the
most simple MSA with sequence identifiers and column names.
"""
abstract type AbstractMultipleSequenceAlignment <: AbstractAlignedObject end

"A MIToS aligned sequence is an `AbstractMatrix{Residue}` with only 1 row/sequence."
abstract type AbstractAlignedSequence <: AbstractAlignedObject end

# Multiple Sequence Alignment
# ===========================

"""
This MSA type include a `NamedArray` wrapping a `Matrix` of `Residue`s. The use of
`NamedArray` allows to store sequence names and original column numbers as `String`s, and
fast indexing using them.
"""
mutable struct MultipleSequenceAlignment <: AbstractMultipleSequenceAlignment
    matrix::NamedArray{ Residue, 2, Array{Residue, 2},
                        Tuple{OrderedDict{String, Int},
                        OrderedDict{String, Int}} }

    function (::Type{MultipleSequenceAlignment})(matrix::NamedArray{Residue,2})
        setdimnames!(matrix,("Seq","Col"))
        new(matrix)
    end
end

"""
This type represent an MSA, similar to `MultipleSequenceAlignment`, but It also stores
`Annotations`. This annotations are used to store residue coordinates (i.e. mapping
to UniProt residue numbers).
"""
mutable struct AnnotatedMultipleSequenceAlignment <: AbstractMultipleSequenceAlignment
    matrix::NamedArray{ Residue, 2, Array{Residue, 2},
                        Tuple{OrderedDict{String, Int},
                        OrderedDict{String, Int}} }
    annotations::Annotations

    function (::Type{AnnotatedMultipleSequenceAlignment})(matrix::NamedArray{Residue,2},
                                                          annotations::Annotations)
        setdimnames!(matrix,("Seq","Col"))
        new(matrix, annotations)
    end
end

# Aligned Sequences
# -----------------

"""
An `AlignedSequence` wraps a `NamedArray{Residue,2}` with only 1 row/sequence. The
`NamedArray` stores the sequence name and original column numbers as `String`s.
"""
mutable struct AlignedSequence <: AbstractAlignedSequence
    matrix::NamedArray{ Residue, 2, Array{Residue, 2},
                        Tuple{OrderedDict{String, Int},
                        OrderedDict{String, Int}} }

    function (::Type{AlignedSequence})(matrix::NamedArray{Residue,2})
        @assert size(matrix,1) == 1 "There are more than one sequence."
        setdimnames!(matrix,("Seq","Col"))
        new(matrix)
    end
end

"""
This type represent an aligned sequence, similar to `AlignedSequence`, but It also stores
its `Annotations`.
"""
mutable struct AnnotatedAlignedSequence <: AbstractAlignedSequence
    matrix::NamedArray{ Residue, 2, Array{Residue, 2},
                        Tuple{OrderedDict{String, Int},
                        OrderedDict{String, Int}} }
    annotations::Annotations

    function (::Type{AnnotatedAlignedSequence})(matrix::NamedArray{Residue,2},
                                                annotations::Annotations)
        @assert size(matrix,1) == 1 "There are more than one sequence."
        setdimnames!(matrix,("Seq","Col"))
        new(matrix, annotations)
    end
end

# Constructors
# ------------

function AnnotatedMultipleSequenceAlignment(msa::NamedArray{Residue,2})
    AnnotatedMultipleSequenceAlignment(msa, Annotations())
end

function AnnotatedMultipleSequenceAlignment(msa::Matrix{Residue})
    AnnotatedMultipleSequenceAlignment(NamedArray(msa))
end

function MultipleSequenceAlignment(msa::Matrix{Residue})
    MultipleSequenceAlignment(NamedArray(msa))
end

function AnnotatedAlignedSequence(seq::NamedArray{Residue,2})
    AnnotatedAlignedSequence(seq, Annotations())
end

function AnnotatedAlignedSequence(seq::Matrix{Residue})
    AnnotatedAlignedSequence(NamedArray(seq))
end

function AlignedSequence(seq::Matrix{Residue})
    AlignedSequence(NamedArray(seq))
end

# AnnotatedAlignedObject
# ----------------------

const AnnotatedAlignedObject = Union{ AnnotatedMultipleSequenceAlignment,
                                        AnnotatedAlignedSequence    }

const UnannotatedAlignedObject = Union{   MultipleSequenceAlignment,
                                            AlignedSequence    }

# Matrices
# --------

const MSAMatrix = Union{ Matrix{Residue}, NamedArray{Residue,2} }

# Getters
# -------

"`annotations` returns the `Annotations` of an MSA or aligned sequence."
@inline annotations(msa::AnnotatedMultipleSequenceAlignment) = msa.annotations
@inline annotations(seq::AnnotatedAlignedSequence) = seq.annotations

"`namedmatrix` returns the `NamedArray{Residue,2}` stored in an MSA or aligned sequence."
@inline namedmatrix(x::AbstractAlignedObject) = x.matrix

# Convert
# -------

function Base.convert(::Type{MultipleSequenceAlignment},
                      msa::AnnotatedMultipleSequenceAlignment)
    MultipleSequenceAlignment(namedmatrix(msa))
end

function Base.convert(::Type{AlignedSequence}, seq::AnnotatedAlignedSequence)
    AlignedSequence(namedmatrix(seq))
end

function Base.convert(::Type{AnnotatedMultipleSequenceAlignment},
                      msa::MultipleSequenceAlignment)
    AnnotatedMultipleSequenceAlignment(namedmatrix(msa), Annotations())
end

function Base.convert(::Type{AnnotatedAlignedSequence}, seq::AlignedSequence)
    AnnotatedAlignedSequence(namedmatrix(seq), Annotations())
end

# AbstractArray Interface
# -----------------------

for f in (:size, :length)
    @eval Base.$(f)(x::AbstractAlignedObject) = $(f)(namedmatrix(x))
end

for T in (  :(AlignedSequence),
            :(AnnotatedAlignedSequence),
            :(MultipleSequenceAlignment),
            :(AnnotatedMultipleSequenceAlignment)  )
    @eval Base.linearindexing(::Type{$(T)}) = Base.LinearFast()
end

@inline Base.getindex(x::AbstractAlignedObject,
                      args...) = getindex(namedmatrix(x), args...)

@inline function Base.setindex!(x::AbstractAlignedObject, value, args...)
    setindex!(namedmatrix(x), value, args...)
end

# Special getindex/setindex! for sequences to avoid `seq["seqname","colname"]`

@inline Base.getindex(x::AbstractAlignedSequence, i) = getindex(namedmatrix(x), 1, i)

@inline function Base.setindex!(x::AbstractAlignedSequence, value, i)
    setindex!(namedmatrix(x), value, 1, i)
end

# Show
# ----

for T in (  :(AlignedSequence),
            :(AnnotatedAlignedSequence),
            :(MultipleSequenceAlignment),
            :(AnnotatedMultipleSequenceAlignment)  )
    @eval begin

        Base.show(io::IO, ::MIME"text/plain", x::$(T)) = show(io, x)

        function Base.show(io::IO, x::$(T))
            type_name = split(string($T),'.')[end]
            if isa(x, AnnotatedAlignedObject)
                print(io, type_name, " with ", length(annotations(x)), " annotations : ")
            else
                print(io, type_name, " : ")
            end
            show(io, namedmatrix(x))
        end

    end
end

# Transpose
# ---------
#
# transpose is ~ 0.00022 seconds faster than ctranspose for PF00085
#

Base.transpose(x::AbstractAlignedObject)  = transpose(namedmatrix(x))
Base.ctranspose(x::AbstractAlignedObject) = transpose(namedmatrix(x))

# Selection without Mappings
# --------------------------

"""
`getresidues` allows you to access the residues stored inside an MSA or aligned sequence
as a `Matrix{Residue}` without annotations nor column/row names.
"""
getresidues(x::Matrix{Residue}) = x
getresidues(x::NamedArray{Residue,2}) = getarray(x)
getresidues(x::AbstractAlignedObject) = getresidues(namedmatrix(x))

"`nsequences` returns the number of sequences on the MSA."
nsequences(x::AbstractMatrix{Residue}) = size(x, 1)

"`ncolumns` returns the number of MSA columns or positions."
ncolumns(x::AbstractMatrix{Residue}) = size(x, 2)

"""
`getresiduesequences` returns a `Vector{Vector{Residue}}` with all the MSA sequences without
annotations nor column/sequence names.
"""
function getresiduesequences(msa::Matrix{Residue})
    nseq = nsequences(msa)
    tmsa = permutedims(msa, [2,1])
    sequences = Array(Vector{Residue}, nseq)
    for i in 1:nseq
        @inbounds sequences[i] = tmsa[:,i]
    end
    sequences
end

getresiduesequences(x::NamedArray{Residue,2}) = getresiduesequences(getresidues(x))
getresiduesequences(x::AbstractAlignedObject) = getresiduesequences(getresidues(x))

# Select sequence
# ---------------

# Gives you the annotations of the Sequence
function getsequence(data::Annotations, id::String)
    GS = Dict{Tuple{String,String},String}()
    GR = Dict{Tuple{String,String},String}()
    if length(data.sequences) > 0
        for (key, value) in data.sequences
            if key[1] == id
                GS[key] = value
            end
        end
        sizehint!(GS, length(GS))
    end
    if length(data.residues) > 0
        for (key, value) in data.residues
            if key[1] == id
                GR[key] = value
            end
        end
        sizehint!(GR, length(GR))
    end
    Annotations(data.file, GS, data.columns, GR)
end

@doc """
`getsequence` takes an MSA and a sequence number or identifier and returns an aligned
sequence object. If the MSA is an `AnnotatedMultipleSequenceAlignment`, it returns an
`AnnotatedAlignedSequence` with the sequence annotations. From a
`MultipleSequenceAlignment`, It returns an `AlignedSequence` object. If an `Annotations`
object and a sequence identifier are used, this function returns the annotations related
to the sequence.
""" getsequence

getsequence(msa::Matrix{Residue}, i::Int) = msa[i:i,:]

getsequence(msa::NamedArray{Residue,2}, i::Int) = msa[i:i,:]
getsequence(msa::NamedArray{Residue,2}, id::String) = msa[String[id],:]

function getsequence(msa::AnnotatedMultipleSequenceAlignment, i::Int)
    seq   = namedmatrix(msa)[i:i,:]
    annot = getsequence(annotations(msa), names(seq, 1)[1])
    AnnotatedAlignedSequence(seq, annot)
end

function getsequence(msa::AnnotatedMultipleSequenceAlignment, id::String)
    seq   = namedmatrix(msa)[String[id],:]
    annot = getsequence(annotations(msa), id)
    AnnotatedAlignedSequence(seq, annot)
end

function getsequence(msa::MultipleSequenceAlignment, seq::Union{Int,String})
    AlignedSequence(getsequence(namedmatrix(msa), seq))
end

# Names
# -----

"""
`sequencenames(msa)`

It returns a `Vector{String}` with the sequence names/identifiers.
"""
sequencenames(x::NamedArray{Residue,2})::Vector{String} = names(x,1)
sequencenames(x::AbstractAlignedObject)::Vector{String} = sequencenames(namedmatrix(x))
sequencenames(msa::AbstractMatrix{Residue})::Vector{String} = map(string, 1:size(msa,1))

"""
`columnnames(msa)`

It returns a `Vector{String}` with the sequence names/identifiers. If the `msa` is a
`Matrix{Residue}` this function returns the actual column numbers as strings. Otherwise it
returns the column number of the original MSA through the wrapped `NamedArray` column names.
"""
columnnames(x::NamedArray{Residue,2})::Vector{String} = names(x,2)
columnnames(x::AbstractAlignedObject)::Vector{String} = columnnames(namedmatrix(x))
columnnames(msa::AbstractMatrix{Residue})::Vector{String} = map(string, 1:size(msa,2))

# Copy, deepcopy
# --------------

for f in (:copy, :deepcopy)
    @eval begin
        function Base.$(f)(msa::AnnotatedMultipleSequenceAlignment)
            AnnotatedMultipleSequenceAlignment( $(f)(namedmatrix(msa)),
                                                $(f)(annotations(msa)) )
        end
        function Base.$(f)(msa::MultipleSequenceAlignment)
            MultipleSequenceAlignment($(f)(namedmatrix(msa)))
        end
        function Base.$(f)(seq::AnnotatedAlignedSequence)
            AnnotatedAlignedSequence($(f)(seq.matrix), $(f)(seq.annotations))
        end
        Base.$(f)(seq::AlignedSequence) = AlignedSequence($(f)(seq.matrix))
    end
end

# Get annotations
# ---------------

for getter in ( :getannotcolumn, :getannotfile, :getannotresidue, :getannotsequence )
    @eval $(getter)(x::AnnotatedAlignedObject, args...) = $(getter)(annotations(x), args...)
end

# Set annotations
# ---------------

for setter in ( :setannotcolumn!, :setannotfile!, :setannotresidue!, :setannotsequence!,
                :annotate_modification!,
                :delete_annotated_modifications!,
                :printmodifications )
    @eval $(setter)(x::AnnotatedAlignedObject, args...) = $(setter)(annotations(x), args...)
end

# To be used on AbstractMultipleSequenceAlignment methods
@inline function annotate_modification!(msa::MultipleSequenceAlignment, str::String)
    # It's generally used in a boolean context: annotate && annotate_modification!(...
    false
end

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

"""
It returns a `Vector{Int}` with the original column number of each column on the actual MSA.
The mapping is annotated in the "ColMap" file annotation of an
`AnnotatedMultipleSequenceAlignment` or in the column names of an `NamedArray` or
`MultipleSequenceAlignment`.
"""
function getcolumnmapping(msa::AnnotatedMultipleSequenceAlignment)
    annot = getannotfile(msa)
    if haskey(annot, "ColMap")
        return _str2int_mapping(annot["ColMap"])
    else
        return getcolumnmapping(namedmatrix(msa))
    end
end

function getcolumnmapping(msa::NamedArray{Residue,2})
    Int[ parse(Int,pos) for pos in names(msa,2) ]
end

getcolumnmapping(msa::MultipleSequenceAlignment) = getcolumnmapping(namedmatrix(msa))

"""
It returns the sequence coordinates as a `Vector{Int}` for an MSA sequence. That vector has
one element for each MSA column. If the number if `0` in the mapping, there is a gap in
that column for that sequence.
"""
function getsequencemapping(msa::AnnotatedMultipleSequenceAlignment, seq_id::String)
    _str2int_mapping(getannotsequence(msa, seq_id, "SeqMap"))
end

function getsequencemapping(msa::AnnotatedMultipleSequenceAlignment, seq_num::Int)
    getsequencemapping(msa, sequencenames(msa)[seq_num])
end

# Sequences as strings
# --------------------

"""
```
stringsequence(seq)
stringsequence(msa, i::Int)
stringsequence(msa, id::String)
```

It returns the selected sequence as a `String`.
"""
stringsequence(msa::AbstractMatrix{Residue}, i) = String(vec(msa[i,:]))

function stringsequence(msa::AbstractMultipleSequenceAlignment, i)
    stringsequence(namedmatrix(msa), i)
end

function stringsequence(seq::AbstractMatrix{Residue})
    @assert size(seq,1) == 1 "There are more than one sequence/row."
    String(vec(seq))
end

stringsequence(seq::AbstractAlignedSequence) = stringsequence(namedmatrix(seq))
