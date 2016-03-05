import Base: print, show, copy, deepcopy, empty!, isempty, print, show

# Annotations
# ===========

"""
The `Annotations` type is basically a container for `Dict`s with the annotations of a multiple sequence alignment.
`Annotations` was designed for storage of annotations of the **Stockholm format**.

MIToS also uses MSA annotations to keep track of:

- **Modifications** of the MSA (`MIToS_...`) as deletion of sequences or columns.
- Positions numbers in the original MSA file (**column mapping:** `ColMap`)
- Position of the residues in the sequence (**sequence mapping:** `SeqMap`)
"""
type Annotations
    file::OrderedDict{ASCIIString, ByteString}
    sequences::Dict{Tuple{ASCIIString,ASCIIString},ASCIIString}
    columns::Dict{ASCIIString,ASCIIString}
    residues::Dict{Tuple{ASCIIString,ASCIIString},ASCIIString}
end

call(::Type{Annotations}) = Annotations(OrderedDict{ASCIIString, ByteString}(),
                                        Dict{Tuple{ASCIIString,ASCIIString},ASCIIString}(),
                                        Dict{ASCIIString, ASCIIString}(),
                                        Dict{Tuple{ASCIIString,ASCIIString},ASCIIString}())

"""
Creates empty MSA `Annotations` of length 0 using `sizehint!`
"""
empty(::Type{Annotations}) = Annotations( sizehint!(OrderedDict{ASCIIString, ByteString}(), 0), # There is not sizehint! for OrderedDict right now
                                         sizehint!(Dict{Tuple{ASCIIString,ASCIIString},ASCIIString}(), 0),
                                         sizehint!(Dict{ASCIIString, ASCIIString}(), 0),
                                         sizehint!(Dict{Tuple{ASCIIString,ASCIIString},ASCIIString}(), 0))

# Filters
# -------

# This function is useful because of the Julia issue 12495
_filter(str::ASCIIString, mask) = ASCIIString( str.data[mask] )

"""
For filter column and sequence mapping of the format: ",,,,10,11,,12"
"""
_filter_mapping(str_map::ASCIIString, mask) = join(split(str_map, ',')[mask], ',')

"""
`filtersequences!(data::Annotations, ids::IndexedArray, mask::AbstractArray{Bool,1})` is useful for deleting annotations for a group of sequences.
`ids` should be an `IndexedArray` with the `seqname`s of the annotated sequences and `mask` should be a logical vector.
"""
function filtersequences!(data::Annotations, ids::IndexedArray, mask::AbstractArray{Bool,1})
    if length(data.sequences) > 0 || length(data.residues) > 0
        del = ids[ !mask ]
        for key in keys(data.residues)
            if key[1] in del
                delete!(data.residues, key)
            end
        end
        for key in keys(data.sequences)
            if key[1] in del
                delete!(data.sequences, key)
            end
        end
        data.sequences = sizehint!(data.sequences, length(data.sequences))
        data.residues = sizehint!(data.residues, length(data.residues))
    end
    data
end

"""
`filtercolumns!(data::Annotations, mask)` is useful for deleting annotations for a group of columns (creating a subset in place).
"""
function filtercolumns!(data::Annotations, mask)
    if length(data.residues) > 0
        for (key,value) in data.residues
            data.residues[key] = _filter(value, mask)
        end
    end
    if length(data.columns) > 0
        for (key,value) in data.columns
            data.columns[key] = _filter(value, mask)
        end
    end
    if length(data.sequences) > 0
        for (key,value) in data.sequences
            if key[2] == "SeqMap"
                data.sequences[key] = _filter_mapping(value, mask)
            end
        end
    end
    filecolumnmapping = get(data.file, "ColMap", "")
    if filecolumnmapping != ""
        data.file["ColMap"] = _filter_mapping(filecolumnmapping, mask)
    end
    data
end

# Copy, deepcopy and empty!
# -------------------------

for fun in [:copy, :deepcopy]
    @eval $(fun)(ann::Annotations) = Annotations( $(fun)( ann.file ), $(fun)( ann.sequences ), $(fun)( ann.columns ), $(fun)( ann.residues ))
end

function empty!(ann::Annotations)
    empty!(ann.file)
    empty!(ann.sequences)
    empty!(ann.columns)
    empty!(ann.residues)
    ann
end

isempty(ann::Annotations) = isempty(ann.file) && isempty(ann.sequences) && isempty(ann.columns) && isempty(ann.residues)

# ncolumns
# --------

"""
`ncolumns(ann::Annotations)` returns the number of columns/residues with annotations.
This function returns `-1` if there is not annotations per column/residue.
"""
function ncolumns(ann::Annotations)
    for value in values(ann.columns)
        return(length(value))
    end
    for value in values(ann.residues)
        return(length(value))
    end
    -1
end

# Getters & Setters
# -----------------

for (fun, field) in [ (:getannotfile,   :(ann.file)),
                     (:getannotcolumn, :(ann.columns))]
    @eval begin
        $(fun)(ann::Annotations) = $(field)
        $(fun)(ann::Annotations, feature::ASCIIString) = getindex($(field), feature)
        $(fun)(ann::Annotations, feature::ASCIIString, default::ASCIIString) = get($(field), feature, default)
    end
end

for (fun, field) in [ (:getannotsequence, :(ann.sequences)),
                     (:getannotresidue,  :(ann.residues))]
    @eval begin
        $(fun)(ann::Annotations) = $(field)
        $(fun)(ann::Annotations, seqname::ASCIIString, feature::ASCIIString) = getindex($(field), (seqname, feature))
        $(fun)(ann::Annotations, seqname::ASCIIString, feature::ASCIIString, default::ASCIIString) = get($(field), (seqname, feature), default)
    end
end

function _test_feature_name(feature::ASCIIString)
    length(feature) <= 50 || throw(ErrorException("Feature name has a limit of 50 characters."))
    ismatch(r"\s", feature) && throw(ErrorException("Feature name must not have spaces."))
end

function setannotfile!(ann::Annotations, feature::ASCIIString, annotation::ByteString)
    _test_feature_name(feature)
    previous = get(ann.file, feature, "")
    ann.file[feature] = previous != "" ? string(previous, '\n', annotation) : annotation
end

function setannotsequence!(ann::Annotations, seqname::ASCIIString, feature::ASCIIString, annotation::ASCIIString)
    _test_feature_name(feature)
    previous = get(ann.sequences, (seqname, feature), "")
    ann.sequences[(seqname, feature)] = previous != "" ? string(previous, '\n', annotation) : annotation
end

function setannotcolumn!(ann::Annotations, feature::ASCIIString, annotation::ASCIIString)
    _test_feature_name(feature)
    len = ncolumns(ann)
    if (len == -1) || (len == length(annotation))
        setindex!(ann.columns, annotation, feature)
    else
        throw(DimensionMismatch("You should have exactly 1 char per column ($len columns/residues)"))
    end
end

function setannotresidue!(ann::Annotations, seqname::ASCIIString, feature::ASCIIString, annotation::ASCIIString)
    _test_feature_name(feature)
    len = ncolumns(ann)
    if (len == -1) || (len == length(annotation))
        setindex!(ann.residues, annotation, (seqname, feature))
    else
        throw(DimensionMismatch("You should have exactly 1 char per residue ($len columns/residues)"))
    end
end

@doc "`getannotfile(ann[, feature[,default]])` returns per file annotation for `feature`" getannotfile
@doc "`getannotcolumn(ann[, feature[,default]])` returns per column annotation for `feature`" getannotcolumn
@doc "`getannotsequence(ann[, seqname, feature[,default]])` returns per sequence annotation for `(seqname, feature)`" getannotsequence
@doc "`getannotresidue(ann[, seqname, feature[,default]])` returns per residue annotation for `(seqname, feature)`" getannotresidue
@doc "`setannotfile!(ann, feature, annotation)` stores per file `annotation` for `feature`" setannotfile!
@doc "`setannotcolumn!(ann, feature, annotation)` stores per column `annotation` (1 char per column) for `feature`" setannotcolumn!
@doc "`setannotsequence!(ann, seqname, feature, annotation)` stores per sequence `annotation` for `(seqname, feature)`" setannotsequence!
@doc "`setannotresidue!(ann, seqname, feature, annotation)` stores per residue `annotation` (1 char per residue) for `(seqname, feature)`" setannotresidue!

# MIToS modification annotations
# ===============================

"""
Annotates on file annotations the modifications realized by MIToS on the MSA
"""
function annotate_modification!(ann::Annotations, modification::ASCIIString)
    setannotfile!(ann, string("MIToS_", Dates.now()), modification)
    true # generally used on bool context: annotate && annotate_modification!(...
end

"""
Deletes all the MIToS annotated modifications
"""
function delete_annotated_modifications!(ann::Annotations)
    for key in keys(ann.file)
        if length(key) > 6 && key[1:6] == "MIToS_"
            delete!(ann.file, key)
        end
    end
end

using Base.Markdown

"""
Prints MIToS annotated modifications
"""
function printmodifications(ann::Annotations)
    for (key,value) in ann.file
        if length(key) > 6 && key[1:6] == "MIToS_"
            list_k = split(key, '_')
            println("-------------------")
            println(list_k[2])
            println('\n', value)
        end
    end
end

# Show & Print Annotations
# ========================

function _printfileannotations(io::IO, ann::Annotations)
    if !isempty(ann.file)
        for (key, value) in ann.file
            for val in split(value, '\n')
                println(io, string("#=GF ", key, "   ", val))
            end
        end
    end
end

function _printcolumnsannotations(io::IO, ann::Annotations)
    if !isempty(ann.columns)
        for (key, value) in ann.columns
            println(io, string("#=GC ", key, "\t\t\t", value))
        end
    end
end

function _printsequencesannotations(io::IO, ann::Annotations)
    if !isempty(ann.sequences)
        for (key, value) in ann.sequences
            for val in split(value, '\n')
                println(io, string("#=GS ", key[1], '\t', key[2], ' ', val))
            end
        end
    end
end

function _printresiduesannotations(io::IO, ann::Annotations)
    if !isempty(ann.residues)
        for (key, value) in ann.residues
            println(io, string("#=GR ", key[1], '\t', key[2], '\t', value))
        end
    end
end

function print(io::IO, ann::Annotations)
    _printfileannotations(io, ann)
    _printsequencesannotations(io, ann)
    _printresiduesannotations(io, ann)
    _printcolumnsannotations(io, ann)
end

show(io::IO, ann::Annotations) = print(io, ann)
