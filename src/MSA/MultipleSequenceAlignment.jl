# import Base: length, getindex, setindex!, size, copy, deepcopy, empty!,
#              convert, transpose, ctranspose, names

"""
MIToS MSA and aligned sequences (aligned objects) are subtypes of `AbstractMatrix{Residue}`,
because MSAs and sequences are stored as `Matrix` of `Residue`s.
"""
abstract AbstractAlignedObject <: AbstractMatrix{Residue}

"""
MSAs are stored as `Matrix{Residue}`. It's possible to use a `NamedArray{Residue,2}` as the
most simple MSA with sequence identifiers and column names.
"""
abstract AbstractMultipleSequenceAlignment <: AbstractAlignedObject

"A MIToS aligned sequence is an `AbstractMatrix{Residue}` with only 1 row/sequence."
abstract AbstractAlignedSequence <: AbstractAlignedObject

# Multiple Sequence Alignment
# ===========================

"""
This MSA type include a `NamedArray` wrapping a `Matrix` of `Residue`s. The use of
`NamedArray` allows to store sequence names and original column numbers as `String`s, and
fast indexing using them.
"""
type MultipleSequenceAlignment <: AbstractMultipleSequenceAlignment
    matrix::NamedArray{ Residue, 2, Array{Residue, 2},
                        Tuple{OrderedDict{String, Int64},
                        OrderedDict{String, Int64}} }
end

"""
This type represent an MSA, similar to `MultipleSequenceAlignment`, but It also stores
`Annotations`. This annotations are used to store residue coordinates (i.e. mapping
to UniProt residue numbers).
"""
type AnnotatedMultipleSequenceAlignment <: AbstractMultipleSequenceAlignment
    matrix::NamedArray{ Residue, 2, Array{Residue, 2},
                        Tuple{OrderedDict{String, Int64},
                        OrderedDict{String, Int64}} }
    annotations::Annotations
end

# Aligned Sequences
# -----------------

"""
An `AlignedSequence` wraps a `NamedArray{Residue,2}` with only 1 row/sequence. The
`NamedArray` stores the sequence name and original column numbers as `String`s.
"""
type AlignedSequence <: AbstractAlignedSequence
    matrix::NamedArray{ Residue, 2, Array{Residue, 2},
                        Tuple{OrderedDict{String, Int64},
                        OrderedDict{String, Int64}} }
end

"""
This type represent an aligned sequence, similar to `AlignedSequence`, but It also stores
its `Annotations`.
"""
type AnnotatedAlignedSequence <: AbstractAlignedSequence
    matrix::NamedArray{ Residue, 2, Array{Residue, 2},
                        Tuple{OrderedDict{String, Int64},
                        OrderedDict{String, Int64}} }
    annotations::Annotations
end

# AnnotatedAlignedObject
# ----------------------

typealias AnnotatedAlignedObject Union{ AnnotatedMultipleSequenceAlignment,
                                        AnnotatedAlignedSequence    }

typealias UnannotatedAlignedObject Union{   AnnotatedMultipleSequenceAlignment,
                                            AnnotatedAlignedSequence    }

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

@inline Base.setindex!(x::AbstractAlignedObject,
                       args...) = setindex!(namedmatrix(x), args...)

# Show
# ----

for T in (  :(AlignedSequence),
            :(AnnotatedAlignedSequence),
            :(MultipleSequenceAlignment),
            :(AnnotatedMultipleSequenceAlignment)  )
    @eval begin

        Base.show(io::IO, ::MIME"text/plain", x::$(T)) = show(io, x)

        function Base.show(io::IO, x::$(T))
            print(io, string($(T)," : "))
            show(io, namedmatrix(x))
        end

    end
end

# Transpose
# ---------
#
# transpose is ~ 0.00022 seconds faster than ctranspose for PF00085
#

Base.transpose(msa::AbstractMultipleSequenceAlignment)  = transpose(namedmatrix(x))
Base.ctranspose(msa::AbstractMultipleSequenceAlignment) = transpose(namedmatrix(x))

# Selection without Mappings
# --------------------------

"""
`getresidues` allows you to access the residues stored inside an MSA or aligned sequence
as a `Matrix{Residue}` without annotations nor column/row names.
"""
getresidues(x::AbstractAlignedObject) = array(namedmatrix(x))

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
    tmsa = msa'
    sequences = Array(Vector{Residue}, nseq)
    for i in 1:nseq
        @inbounds sequences[i] = tmsa[:,i]
    end
    sequences
end

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

getsequence(msa::AbstractMatrix{Residue}, i::Int) = msa[i:i,:]

function getsequence(msa::AnnotatedMultipleSequenceAlignment, i::Int)
    seq   = namedmatrix(msa)[i:i,:]
    annot = getsequence(annotations(msa), names(seq, 1)[1])
    AnnotatedAlignedSequence(seq, annot)
end

function getsequence(msa::MultipleSequenceAlignment, i::Int)
    AlignedSequence(namedmatrix(msa)[i:i,:])
end

function getsequence(msa::AnnotatedMultipleSequenceAlignment, id::String)
    seq   = namedmatrix(msa)[String[id],:]
    annot = getsequence(annotations(msa), id)
    AnnotatedAlignedSequence(seq, annot)
end

function getsequence(msa::MultipleSequenceAlignment, id::String)
    AlignedSequence(namedmatrix(msa)[String[id],:])
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
