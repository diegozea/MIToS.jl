@doc """
    AnnotatedFASTASequences

Sequence format for:

  - Annotated FASTA Format (AFF)
  - CAID modified FASTA references/labels

References:

  - AFF README (pinned commit):
    `https://github.com/NawarMalhis/AFF/blob/f0e7033bc4d08b255cd883c99523d9379a9fbb16/README.md`
  - CAID modified FASTA description:
    `https://caid.idpcentral.org/challenge/results`

Parser behavior:

  - AFF header lines (`# ...`) before the first sequence are optional.
  - Header key/value metadata is parsed using `key: value` and stored as file
    annotations.
  - The AFF `Format` block defines residue annotation feature names by order.
  - CAID-like files (no header, one annotation line per sequence) store the
    residue annotation as feature `"Feature"`.
  - If there is no header and multiple annotation lines are present, default
    names are `"Feature_1"`, `"Feature_2"`, etc.
  - Keyword `format_has_sequence_line` controls whether the `Format` section
    includes a sequence descriptor line after the `>...` prototype and before
    annotation tags (`true` by default).
"""
struct AnnotatedFASTASequences <: SequenceFormat end

# Header parsing
# ==============

# AFF README uses `key: value` for header metadata.
function _parse_aff_key_value(content::AbstractString)
    parts = split(content, ':'; limit = 2)
    length(parts) == 2 || return nothing
    key = strip(parts[1])
    value = strip(parts[2])
    isempty(key) && return nothing
    return key, value
end

"""
    _flush_aff_multiline_section!(
        file_annotations,
        current_multiline_key,
        current_multiline_value,
    ) -> AbstractString

If a multiline AFF header section is open, append its collected lines into
`file_annotations[current_multiline_key]` using newline separators and clear
`current_multiline_value`.

Returns the next section key state:

  - `""` when a section was flushed
  - `current_multiline_key` unchanged when there is no open section
"""
function _flush_aff_multiline_section!(
    file_annotations::OrderedDict{String,String},
    current_multiline_key::AbstractString,
    current_multiline_value::Vector{String},
)
    if !isempty(current_multiline_key)
        _append_with_separator!(
            file_annotations,
            current_multiline_key,
            join(current_multiline_value, "\n"),
            "\n",
        )
        empty!(current_multiline_value)
        return ""
    end
    return current_multiline_key
end

function _parse_aff_header(
    header_lines::Vector{String},
    ;
    format_has_sequence_line::Bool = true,
)::Tuple{OrderedDict{String,String},Vector{String}}
    file_annotations = OrderedDict{String,String}()
    format_tags = String[]
    current_multiline_key = ""
    current_multiline_value = String[]
    comment_counter = 0
    # AFF header has a single Format section.
    format_line_counter = 0

    for line in header_lines
        startswith(line, '#') || continue
        content = strip(line[2:end])
        keyvalue = isempty(content) ? nothing : _parse_aff_key_value(content)
        section_boundary = isempty(content) || keyvalue !== nothing
        if section_boundary
            # A multiline section can end with either an empty `#` line or a new key.
            current_multiline_key = _flush_aff_multiline_section!(
                file_annotations,
                current_multiline_key,
                current_multiline_value,
            )
        end

        isempty(content) && continue

        if keyvalue !== nothing
            key, value = keyvalue

            if isempty(value)
                current_multiline_key = String(key)
            else
                _append_with_separator!(file_annotations, key, value, "\n")
            end
            continue
        end

        # Content lines that belong to a multiline section.
        if !isempty(current_multiline_key)
            push!(current_multiline_value, content)

            # The `Format` section lines define annotation feature names.
            if lowercase(current_multiline_key) == "format"
                format_line_counter += 1
                # Create `format_tags` to store the feature names.
                # For `format_tags`, skip:
                #   1. the first line in Format (`>...` prototype)
                if format_line_counter == 1
                    continue
                end
                #   2. the second line when `format_has_sequence_line` is true
                if format_has_sequence_line && format_line_counter == 2
                    continue
                end
                #   3. any empty lines
                tag = split(content, r"\s+"; limit = 2)[1]
                isempty(tag) || push!(format_tags, tag)
            end
        else
            # Optional top/bottom comment lines outside named sections.
            comment_counter += 1
            file_annotations["Comment_$(comment_counter)"] = content
        end
    end

    # Flush final multiline section.
    current_multiline_key = _flush_aff_multiline_section!(
        file_annotations,
        current_multiline_key,
        current_multiline_value,
    )

    length(unique(format_tags)) == length(format_tags) ||
        throw(ErrorException("The AFF header Format section has repeated feature names."))

    return file_annotations, format_tags
end

# Record parsing
# ==============

function _resolve_annotated_fasta_feature_names(
    n_annotations::Int,
    format_tags::Vector{String},
    sequence_id::String,
)
    if !isempty(format_tags)
        if length(format_tags) != n_annotations
            throw(
                ErrorException(
                    "The sequence $sequence_id has $n_annotations annotation lines, but " *
                    "the header Format section defines $(length(format_tags)) features.",
                ),
            )
        end
        return format_tags
    end

    if n_annotations == 1
        return ["Feature"]
    end

    return ["Feature_$(i)" for i = 1:n_annotations]
end

function _parse_annotated_fasta_record_lines(
    lines::Vector{String},
    sequence_id::String,
    format_tags::Vector{String},
)
    # remove end-of-line markers, keeping spaces as meaningful characters as they can 
    # possibly be used for residue-level annotations (including trailing spaces)
    map!(chomp, lines, lines)
    # remove truly empty lines
    filter!(line -> !isempty(line), lines)

    if isempty(lines)
        throw(ErrorException("There are no sequence lines for $sequence_id."))
    elseif length(lines) == 1
        throw(ErrorException("There are no annotation lines for $sequence_id."))
    end

    # AFF record structure: first line is the sequence, then annotations
    sequence = lines[1]
    annotation_lines = @view lines[2:end]

    # all annotation strings must align residue-by-residue to the sequence
    expected_length = length(sequence)
    for (i, annotation) in enumerate(annotation_lines)
        if length(annotation) != expected_length
            throw(
                ErrorException(
                    "The annotation line $(i) for $sequence_id has " *
                    "$(length(annotation)) characters. $expected_length are expected.",
                ),
            )
        end
    end

    # resolve annotation feature names from AFF Format tags or use defaults
    feature_names = _resolve_annotated_fasta_feature_names(
        length(annotation_lines),
        format_tags,
        sequence_id,
    )

    # Each feature maps to exactly one annotation line in AFF/CAID records.
    residue_annotations = OrderedDict{String,String}()
    for (feature, annotation) in zip(feature_names, annotation_lines)
        residue_annotations[feature] = annotation
    end

    return sequence, residue_annotations
end

function _push_annotated_fasta_record!(
    ids::Vector{String},
    seqs::Vector{String},
    gs::Dict{Tuple{String,String},String},
    gr::Dict{Tuple{String,String},String},
    disambiguator::OnlineSequenceNameDisambiguator,
    original_id::AbstractString,
    lines::Vector{String},
    format_tags::Vector{String},
)
    original_name = String(original_id)
    new_name = _disambiguate_seqname!(disambiguator, original_name)
    push!(ids, new_name)

    if new_name != original_name
        gs[(new_name, "OriginalSeqName")] = original_name
    end

    sequence, residue_annotations =
        _parse_annotated_fasta_record_lines(lines, new_name, format_tags)
    push!(seqs, sequence)

    for (feature, annotation) in residue_annotations
        gr[(new_name, feature)] = annotation
    end

    return nothing
end

# File parser
# ===========

function _pre_read_annotated_fasta_sequences(
    io::Union{IO,AbstractString};
    format_has_sequence_line::Bool = true,
)
    ids = String[]
    seqs = String[]
    gf = OrderedDict{String,String}()
    gs = Dict{Tuple{String,String},String}()
    gr = Dict{Tuple{String,String},String}()
    disambiguator = OnlineSequenceNameDisambiguator()

    header_lines = String[]
    format_tags = String[]
    parsed_header = false

    current_id = ""
    current_lines = String[]

    for line::String in lineiterator(io)
        line = chomp(line)
        isempty(line) && continue

        # store header lines until the first sequence is found defining a current_id
        if startswith(line, '#')
            if !parsed_header && isempty(current_id)
                push!(header_lines, line)
            end
            continue
        end

        if startswith(line, '>') # new sequence record starts, current_id will be updated

            # we parse the header as soon as we find the first sequence
            if !parsed_header
                gf, format_tags = _parse_aff_header(
                    header_lines;
                    format_has_sequence_line = format_has_sequence_line,
                )
                parsed_header = true
            end

            # before updating current_id, we push the previous record if it exists
            if !isempty(current_id)
                _push_annotated_fasta_record!(
                    ids,
                    seqs,
                    gs,
                    gr,
                    disambiguator,
                    current_id,
                    current_lines,
                    format_tags,
                )
            end

            current_id = String(strip(line[2:end]))
            current_lines = empty!(current_lines)
            continue
        end

        isempty(current_id) || push!(current_lines, line)
    end

    # after the loop, we need to push the last record if it exists
    if !isempty(current_id)
        _push_annotated_fasta_record!(
            ids,
            seqs,
            gs,
            gr,
            disambiguator,
            current_id,
            current_lines,
            format_tags,
        )
    end

    return ids, seqs, gf, gs, gr
end

# Loader
# ======

function _load_sequences(
    io::Union{IO,AbstractString},
    format::Type{AnnotatedFASTASequences};
    create_annotations::Bool = false,
    format_has_sequence_line::Bool = true,
)
    ids, seqs, gf, gs, gr = _pre_read_annotated_fasta_sequences(
        io;
        format_has_sequence_line = format_has_sequence_line,
    )
    annot = Annotations()
    if create_annotations
        annot.file = gf
        annot.sequences = gs
        annot.residues = gr
    end
    return ids, seqs, annot
end

# Parse File
# ==========

function Utils.parse_file(
    io::Union{IO,AbstractString},
    format::Type{AnnotatedFASTASequences};
    generatemapping::Bool = false,
    useidcoordinates::Bool = false,
    deletefullgaps::Bool = true,
    keepinserts::Bool = false,
    format_has_sequence_line::Bool = true,
)::Vector{AnnotatedSequence}
    ids, seqs, annot = _load_sequences(
        io,
        format;
        create_annotations = true,
        format_has_sequence_line = format_has_sequence_line,
    )
    _generate_sequences(ids, seqs, annot)
end
