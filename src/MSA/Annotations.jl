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
    file::OrderedDict{String,String}
    sequences::Dict{Tuple{String,String},String}
    columns::Dict{String,String}
    residues::Dict{Tuple{String,String},String}
end

"""
    _rename_columns(
        annotations::Annotations,
        oldnames::Vector{T},
        newnames::Vector{T},
    ) where {T<:AbstractString}

Update column-related annotations in `annotations` when column names change.

Changed columns are recorded in a file annotation named `"N_ColChanges"` where
`N` is the number of times columns were renamed. Each value is a comma-separated
list of `"old=>new"` pairs, where both names appear quoted as Julia strings. The modification is also tracked with
`annotate_modification!`.
"""
function _rename_columns(
    annotations::Annotations,
    oldnames::Vector{T},
    newnames::Vector{T},
) where {T<:AbstractString}
    new_annotations = copy(annotations)
    if oldnames != newnames
        changes =
            [string(repr(o), "=>", repr(n)) for (o, n) in zip(oldnames, newnames) if o != n]
        if !isempty(changes)
            regex = r"^(\d+)_ColChanges$"
            nums = Int[
                parse(Int, m.captures[1]) for k in keys(new_annotations.file) for
                m in eachmatch(regex, k)
            ]
            next = isempty(nums) ? 1 : maximum(nums) + 1
            setannotfile!(new_annotations, string(next, "_ColChanges"), join(changes, ','))
            annotate_modification!(new_annotations, "rename_columns!")
        end
    end
    new_annotations
end

Annotations() = Annotations(
    OrderedDict{String,String}(),
    Dict{Tuple{String,String},String}(),
    Dict{String,String}(),
    Dict{Tuple{String,String},String}(),
)

# Length
# ------

function Base.length(a::Annotations)
    length(a.file) + length(a.sequences) + length(a.columns) + length(a.residues)
end

# Filters
# -------

# This function is useful because of the Julia issue #12495
function _filter(str::AbstractString, mask::AbstractArray{Bool})
    @assert length(str) == length(mask) "The string and the mask must have the same length"
    buffer = IOBuffer(; sizehint = ncodeunits(str))
    i  = 0
    @inbounds for c in str
        i += 1
        if mask[i]
            write(buffer, c)
        end
    end
    String(take!(buffer))
end

function _filter(str::AbstractString, indexes::AbstractArray{Int})
    isempty(indexes) && return ""
    if isascii(str)
        str[indexes] # fast path: ASCII => indices are valid string indices
    else
        chars = collect(str) # one Vector{Char}
        io = IOBuffer(; sizehint = length(indexes))
        @inbounds for i in indexes
            write(io, chars[i])
        end
        String(take!(io))
    end
end

function _filter_mapping(str_map::AbstractString, mask::AbstractVector{Bool})
    io = IOBuffer(; sizehint = ncodeunits(str_map))
    wrote = false
    nfields = 0
    for part in eachsplit(str_map, ','; keepempty = true)
        nfields += 1
        if nfields > length(mask)
            throw(DimensionMismatch("mask length $(length(mask)) < $nfields fields"))
        end
        if @inbounds mask[nfields]
            wrote && write(io, ',')
            write(io, part)
            wrote = true
        end
    end
    if nfields != length(mask)
        throw(DimensionMismatch("mask length $(length(mask)) > $nfields fields"))
    end
    String(take!(io))
end

# Compute start/stop indices for all comma-separated fields in `str_map`.
function _field_ranges(str_map::AbstractString)
    starts = Int[]
    stops = Int[]

    last = lastindex(str_map)
    start = firstindex(str_map)
    i = start

    while i <= last
        if @inbounds str_map[i] == ','
            push!(starts, start)
            push!(stops, start < i ? prevind(str_map, i) : start - 1)
            start = nextind(str_map, i)
        end
        i = nextind(str_map, i)
    end

    # final field (possibly empty)
    push!(starts, start)
    push!(stops, start <= last ? last : start - 1)

    return starts, stops
end

function _filter_mapping(str_map::AbstractString, mask::AbstractVector{<:Integer})
    isempty(mask) && return ""

    starts, stops = _field_ranges(str_map)
    nentries = length(starts)

    io = IOBuffer(; sizehint = ncodeunits(str_map))
    first = true

    for raw_idx in mask
        idx = Int(raw_idx)
        1 <= idx <= nentries || throw(BoundsError(mask, raw_idx))

        first || write(io, ',')
        first = false

        s = @inbounds starts[idx]
        e = @inbounds stops[idx]
        s <= e && write(io, SubString(str_map, s, e))  # empty field ⇒ nothing
    end

    return String(take!(io))
end

"""
`filtersequences!(data::Annotations, ids::Vector{String}, mask::AbstractArray{Bool,1})`

It is useful for deleting sequence annotations. `ids` should be a list of the sequence
names and `mask` should be a logical vector.
"""
function filtersequences!(
    data::Annotations,
    ids::Vector{String},
    mask::AbstractVector{Bool},
)
    @assert length(ids) == length(mask) "It's needed one sequence id per element in the mask."
    nresannot = length(data.residues)
    nseqannot = length(data.sequences)
    if nresannot > 0 || nseqannot > 0
        del = Set(ids[.!mask])
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

function _get_selected_indexes(input_mask)
    if eltype(input_mask) <: Bool
        return findall(input_mask)
    end
    input_mask
end

"""
`filtercolumns!(data::Annotations, mask)`

It is useful for deleting column annotations (creating a subset in place).
"""
function filtercolumns!(data::Annotations, input_mask)
    mask = _get_selected_indexes(input_mask)
    if length(data.residues) > 0
        for (key, value) in data.residues
            data.residues[key] = _filter(value, mask)
        end
    end
    if length(data.columns) > 0
        for (key, value) in data.columns
            data.columns[key] = _filter(value, mask)
        end
    end
    if length(data.sequences) > 0
        for (key, value) in data.sequences
            if key[2] == "SeqMap"
                data.sequences[key] = _filter_mapping(value, mask)
            end
        end
    end
    # vcat can create multiple ColMap annotations in a single MSA; update all of them
    colmap_keys = filter(endswith("ColMap"), keys(data.file))
    for key in colmap_keys
        data.file[key] = _filter_mapping(data.file[key], mask)
    end
    data
end

# Copy, deepcopy and empty!
# -------------------------

for fun in [:copy, :deepcopy]
    @eval begin
        Base.$(fun)(ann::Annotations) = Annotations(
            $(fun)(ann.file),
            $(fun)(ann.sequences),
            $(fun)(ann.columns),
            $(fun)(ann.residues),
        )
    end
end

function Base.empty!(ann::Annotations)
    empty!(ann.file)
    empty!(ann.sequences)
    empty!(ann.columns)
    empty!(ann.residues)
    ann
end

Base.isempty(ann::Annotations) =
    isempty(ann.file) &&
    isempty(ann.sequences) &&
    isempty(ann.columns) &&
    isempty(ann.residues)

# merge! and merge
# ----------------

const _MERGE_NOTE = md"""
NOTE: This function does not check for consistency among the various `Annotations`.
For example, it does not verify that `SeqMap` annotations have consistent lengths.
"""

"""
    merge!(target::Annotations, sources::Annotations...)

Merge one or more source `Annotations` into a target `Annotations`. This function calls
`Base.merge!` on each of the four dictionaries in the `Annotations` type: `file`,
`sequences`, `columns`, and `residues`. Consequently, it behaves like `Base.merge!` for
dictionaries; if the same key exists in different `Annotations` objects, the value from the
last one is used.

$_MERGE_NOTE
"""
function Base.merge!(target::Annotations, sources::Annotations...)
    for source in sources
        merge!(target.file, source.file)
        merge!(target.sequences, source.sequences)
        merge!(target.columns, source.columns)
        merge!(target.residues, source.residues)
    end

    target
end

"""
    merge(target::Annotations, sources::Annotations...)

Create a new `Annotations` object by merging two or more `Annotations`. If the same
annotation exists in different `Annotations` objects, the value from the last one is used.
See also `merge!`.

$_MERGE_NOTE
"""
Base.merge(target::Annotations, sources::Annotations...) =
    merge!(deepcopy(target), sources...)

# ncolumns
# --------

"""
`ncolumns(ann::Annotations)` returns the number of columns/residues with annotations.
This function returns `-1` if there is not annotations per column/residue.
"""
function ncolumns(ann::Annotations)
    for value in values(ann.columns)
        return (length(value))
    end
    for value in values(ann.residues)
        return (length(value))
    end
    -1
end

# Getters
# -------

for (fun, field) in [(:getannotfile, :(ann.file)), (:getannotcolumn, :(ann.columns))]
    @eval begin
        $(fun)(ann::Annotations) = $(field)
        $(fun)(ann::Annotations, feature::String) = getindex($(field), feature)
        $(fun)(ann::Annotations, feature::String, default::String) =
            get($(field), feature, default)
    end
end

for (fun, field) in
    [(:getannotsequence, :(ann.sequences)), (:getannotresidue, :(ann.residues))]
    @eval begin
        $(fun)(ann::Annotations) = $(field)
        $(fun)(ann::Annotations, seqname::String, feature::String) =
            getindex($(field), (seqname, feature))
        $(fun)(ann::Annotations, seqname::String, feature::String, default::String) =
            get($(field), (seqname, feature), default)
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

function setannotsequence!(
    ann::Annotations,
    seqname::String,
    feature::String,
    annotation::String,
)
    _test_feature_name(feature)
    previous = get(ann.sequences, (seqname, feature), "")
    ann.sequences[(seqname, feature)] =
        previous != "" ? string(previous, '\n', annotation) : annotation
end

function setannotcolumn!(ann::Annotations, feature::String, annotation::String)
    _test_feature_name(feature)
    len = ncolumns(ann)
    if (len == -1) || (len == length(annotation))
        setindex!(ann.columns, annotation, feature)
    else
        throw(
            DimensionMismatch(
                string(
                    "You should have exactly 1 char per column (",
                    len,
                    " columns/residues)",
                ),
            ),
        )
    end
end

function setannotresidue!(
    ann::Annotations,
    seqname::String,
    feature::String,
    annotation::String,
)
    _test_feature_name(feature)
    len = ncolumns(ann)
    if (len == -1) || (len == length(annotation))
        setindex!(ann.residues, annotation, (seqname, feature))
    else
        throw(
            DimensionMismatch(
                string(
                    "You should have exactly 1 char per residue (",
                    len,
                    " columns/residues)",
                ),
            ),
        )
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

"""
Deletes all the MIToS annotated modifications
"""
function delete_annotated_modifications!(ann::Annotations)
    for key in keys(ann.file)
        if startswith(key, "MIToS_")
            delete!(ann.file, key)
        end
    end
end

"""
Prints MIToS annotated modifications
"""
function printmodifications(io::IO, ann::Annotations)
    for (key, value) in ann.file
        if startswith(key, "MIToS_")
            list_k = split(key, '_')
            println(io, "-----------------------")
            println(io, list_k[2])
            println(io)
            println(io, value)
        end
    end
end

printmodifications(ann::Annotations) = printmodifications(stdout, ann)

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

# Rename sequences
# ================

# This private function is used by rename sequences during MSA concatenation,
# but it could be useful for other purposes

"""
    _rename_sequences(annotations::Annotations, old2new::Dict{String,String})

Renames sequences in a given `Annotations` object based on a mapping provided in `old2new`.

This function iterates through residue-level and sequence-level annotations in the
provided `Annotations` object. For each annotation, it checks if the sequence name exists
in the `old2new` dictionary. If so, the sequence name is updated to the new name from
`old2new`; otherwise, the original sequence name is retained.

The function then returns a new `Annotations` object with updated sequence names.
"""
function _rename_sequences(annotations::Annotations, old2new::Dict{String,String})
    res_annotations = Dict{Tuple{String,String},String}()
    for ((seqname, annot_name), value) in getannotresidue(annotations)
        seqname = get(old2new, seqname, seqname)
        res_annotations[(seqname, annot_name)] = value
    end
    seq_annotations = Dict{Tuple{String,String},String}()
    for ((seqname, annot_name), value) in getannotsequence(annotations)
        seqname = get(old2new, seqname, seqname)
        seq_annotations[(seqname, annot_name)] = value
    end
    Annotations(
        copy(getannotfile(annotations)),
        seq_annotations,
        copy(getannotcolumn(annotations)),
        res_annotations,
    )
end
