# Annotations
# ===========

"""
The `Annotations` type is basically a container for `Dict`s with the annotations of a
multiple sequence alignment. `Annotations` was designed for storage of annotations of
the **Stockholm format**.

MIToS also uses MSA annotations to keep track of:

- **Modifications** of the MSA (`MIToS_...`) as deletion of sequences or columns.
- Positions numbers in the original MSA file (**column mapping:** `ColMap`)
- Position of the residues in the sequence (**sequence mapping:** `SeqMap`)
"""
@auto_hash_equals mutable struct Annotations
    file::OrderedDict{String, String}
    sequences::Dict{Tuple{String,String},String}
    columns::Dict{String,String}
    residues::Dict{Tuple{String,String},String}
end

Annotations() = Annotations(  OrderedDict{String, String}(),
                              Dict{Tuple{String,String},String}(),
                              Dict{String, String}(),
                              Dict{Tuple{String,String},String}() )

# Length
# ------

function Base.length(a::Annotations)
    length(a.file) + length(a.sequences) + length(a.columns) + length(a.residues)
end

# Filters
# -------

# This function is useful because of the Julia issue #12495
function _filter(str::String, mask::AbstractArray{Bool})
    @assert length(str) == length(mask) "The string and the mask must have the same length"
    #                 data                             readable   writable
    buffer = IOBuffer(Array{UInt8}(undef, lastindex(str)), read=true, write=true)
    # To start at the beginning of the buffer:
    truncate(buffer,0)
    i = 1
    for char in str
        @inbounds if mask[i]
            write(buffer, char)
        end
        i += 1
    end
    String(take!(buffer))
end

_filter(str::String, indexes::AbstractArray{Int}) = String(collect(str)[indexes])

"""
For filter column and sequence mapping of the format: ",,,,10,11,,12"
"""
_filter_mapping(str_map::String, mask) = join(split(str_map, ',')[mask], ',')

"""
`filtersequences!(data::Annotations, ids::Vector{String}, mask::AbstractArray{Bool,1})`

It is useful for deleting sequence annotations. `ids` should be a list of the sequence
names and `mask` should be a logical vector.
"""
function filtersequences!(data::Annotations, ids::Vector{String},
                          mask::AbstractVector{Bool})
    @assert length(ids)==length(mask) "It's needed one sequence id per element in the mask."
    nresannot = length(data.residues)
    nseqannot = length(data.residues)
    if nresannot > 0 || nseqannot > 0
        del = ids[ .!mask ]
    end
    if nresannot > 0
        for key in keys(data.residues)
            if key[1] in del
                delete!(data.residues, key)
            end
        end
        data.residues = sizehint!(data.residues, length(data.residues))
    end
    if nseqannot > 0
        for key in keys(data.sequences)
            if key[1] in del
                delete!(data.sequences, key)
            end
        end
        data.sequences = sizehint!(data.sequences, length(data.sequences))
    end
    data
end

"""
`filtercolumns!(data::Annotations, mask)`

It is useful for deleting column annotations (creating a subset in place).
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
    @eval begin
        Base.$(fun)(ann::Annotations) = Annotations( $(fun)( ann.file ),
                                                     $(fun)( ann.sequences ),
                                                     $(fun)( ann.columns ),
                                                     $(fun)( ann.residues ) )
    end
end

function Base.empty!(ann::Annotations)
    empty!(ann.file)
    empty!(ann.sequences)
    empty!(ann.columns)
    empty!(ann.residues)
    ann
end

Base.isempty(ann::Annotations) = isempty(ann.file) && isempty(ann.sequences) &&
                                 isempty(ann.columns) && isempty(ann.residues)

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

# Getters
# -------

for (fun, field) in [ (:getannotfile,   :(ann.file)),
                      (:getannotcolumn, :(ann.columns))]
    @eval begin
        $(fun)(ann::Annotations) = $(field)
        $(fun)(ann::Annotations, feature::String) = getindex($(field), feature)
        $(fun)(ann::Annotations, feature::String,
               default::String) = get($(field), feature, default)
    end
end

for (fun, field) in [ (:getannotsequence, :(ann.sequences)),
                      (:getannotresidue,  :(ann.residues))]
    @eval begin
        $(fun)(ann::Annotations) = $(field)
        $(fun)(ann::Annotations, seqname::String,
               feature::String) = getindex($(field), (seqname,feature))
        $(fun)(ann::Annotations, seqname::String,
               feature::String, default::String) = get($(field), (seqname,feature), default)
    end
end

@doc """`getannotfile(ann[, feature[,default]])`

It returns per file annotation for `feature`
""" getannotfile

@doc """`getannotcolumn(ann[, feature[,default]])`

It returns per column annotation for `feature`
""" getannotcolumn

@doc """`getannotsequence(ann[, seqname, feature[,default]])`

It returns per sequence annotation for `(seqname, feature)`
""" getannotsequence

@doc """`getannotresidue(ann[, seqname, feature[,default]])`

It returns per residue annotation for `(seqname, feature)`
""" getannotresidue

# Setters
# -------

function _test_feature_name(feature::String)
    @assert length(feature) <= 50 "Feature name has a limit of 50 characters."
    @assert !occursin(r"\s", feature) "Feature name must not have spaces."
end

function setannotfile!(ann::Annotations, feature::String, annotation::String)
    _test_feature_name(feature)
    previous = get(ann.file, feature, "")
    ann.file[feature] = previous != "" ? string(previous, '\n', annotation) : annotation
end

function setannotsequence!(ann::Annotations, seqname::String, feature::String,
                           annotation::String)
    _test_feature_name(feature)
    previous = get(ann.sequences, (seqname, feature), "")
    ann.sequences[(seqname, feature)] = previous != "" ?
                                        string(previous, '\n', annotation) : annotation
end

function setannotcolumn!(ann::Annotations, feature::String, annotation::String)
    _test_feature_name(feature)
    len = ncolumns(ann)
    if (len == -1) || (len == length(annotation))
        setindex!(ann.columns, annotation, feature)
    else
        throw(DimensionMismatch(string("You should have exactly 1 char per column (",
        len, " columns/residues)")))
    end
end

function setannotresidue!(ann::Annotations, seqname::String, feature::String,
                          annotation::String)
    _test_feature_name(feature)
    len = ncolumns(ann)
    if (len == -1) || (len == length(annotation))
        setindex!(ann.residues, annotation, (seqname, feature))
    else
        throw(DimensionMismatch(string("You should have exactly 1 char per residue (",
        len, " columns/residues)")))
    end
end

@doc """`setannotfile!(ann, feature, annotation)`

It stores per file `annotation` for `feature`
""" setannotfile!

@doc """`setannotcolumn!(ann, feature, annotation)`

It stores per column `annotation` (1 char per column) for `feature`
""" setannotcolumn!

@doc """`setannotsequence!(ann, seqname, feature, annotation)`

It stores per sequence `annotation` for `(seqname, feature)`
""" setannotsequence!

@doc """`setannotresidue!(ann, seqname, feature, annotation)`

It stores per residue `annotation` (1 char per residue) for `(seqname, feature)`
""" setannotresidue!

# MIToS modification annotations
# ===============================

"""
Annotates on file annotations the modifications realized by MIToS on the MSA. It always
returns `true`, so It can be used in a boolean context.
"""
function annotate_modification!(ann::Annotations, modification::String)
    setannotfile!(ann, string("MIToS_", Dates.now()), modification)
    true # generally used in a boolean context: annotate && annotate_modification!(...
end

"Deletes all the MIToS annotated modifications"
function delete_annotated_modifications!(ann::Annotations)
    for key in keys(ann.file)
        if startswith(key, "MIToS_")
            delete!(ann.file, key)
        end
    end
end

"Prints MIToS annotated modifications"
function printmodifications(ann::Annotations)
    for (key,value) in ann.file
        if startswith(key, "MIToS_")
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
                println(io, string("#=GF ", key, '\t', val))
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
                println(io, string("#=GS ", key[1], '\t', key[2], '\t', val))
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

function Base.print(io::IO, ann::Annotations)
    _printfileannotations(io, ann)
    _printsequencesannotations(io, ann)
    _printresiduesannotations(io, ann)
    _printcolumnsannotations(io, ann)
end

Base.show(io::IO, ann::Annotations) = print(io, ann)
