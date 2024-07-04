# Sequence File Formats
# =====================

struct FASTASequences <: SequenceFormat end
struct PIRSequences <: SequenceFormat end
struct RawSequences <: SequenceFormat end

_format_fallback(::Type{FASTASequences}) = FASTA
_format_fallback(::Type{PIRSequences}) = PIR
_format_fallback(::Type{RawSequences}) = Raw

function _generate_sequences(ids, seqs, annot)
    n = length(ids)
    @assert n == length(seqs)
    output = AnnotatedSequence[]
    sizehint!(output, n)
    for (id, seq) in zip(ids, seqs)
        seq_annot = getsequence(annot, id)
        push!(output, AnnotatedSequence(id, seq, seq_annot))
    end
    output
end

# Parse File
# ==========

function Utils.parse_file(
    io::Union{IO,AbstractString},
    format::Type{T};
    generatemapping::Bool = false,
    useidcoordinates::Bool = false,
    deletefullgaps::Bool = true,
    keepinserts::Bool = false,
) where {T<:SequenceFormat}
    pre_parser_format = _format_fallback(T)
    ids, seqs, annot = _load_sequences(io, pre_parser_format; create_annotations = true)
    _generate_sequences(ids, seqs, annot)
end

# Print File
# ==========

# PIRSequences
# ------------

# It uses _get_pir_annotations and _print_pir_seq from src/MSA/PIR.jl
function Utils.print_file(
    io::IO,
    seqs::Vector{AnnotatedSequence},
    format::Type{PIRSequences},
)
    for seq in seqs
        seqann = getannotsequence(seq)
        seq_id = sequence_id(seq)
        seq_type, seq_title = _get_pir_annotations(seqann, seq_id)
        seq_str = join(seq)
        _print_pir_seq(io, seq_type, seq_id, seq_title, seq_str)
    end
end

# RawSequences
# ------------

function Utils.print_file(
    io::IO,
    seqs::Vector{AnnotatedSequence},
    format::Type{RawSequences},
)
    for seq in seqs
        println(io, join(seq))
    end
end

# FASTASequences
# --------------

function Utils.print_file(
    io::IO,
    seqs::Vector{AnnotatedSequence},
    format::Type{FASTASequences},
)
    for seq in seqs
        println(io, ">", sequence_id(seq))
        println(io, join(seq))
    end
end
