# Sequence File Formats
# =====================

"""
    FASTASequences <: SequenceFormat

A format type for reading a collection of sequences in FASTA format.
Each sequence defined in the FASTA file is parsed as an individual `AnnotatedSequence`.
This contrasts with the `FASTA` type (which is an `MSAFormat`), where a FASTA file
is interpreted as a single multiple sequence alignment.

Use `FASTASequences` when you want to process a list of sequences from a FASTA file,
rather than treating them as an aligned MSA.
"""
struct FASTASequences <: SequenceFormat end

"""
    PIRSequences <: SequenceFormat

A format type for reading a collection of sequences in PIR/NBRF format.
Each sequence entry in the PIR file is parsed as an individual `AnnotatedSequence`.
This contrasts with the `PIR` type (which is an `MSAFormat`), where a PIR file
is interpreted as a single multiple sequence alignment.

Use `PIRSequences` when you want to process a list of sequences from a PIR file,
rather than treating them as an aligned MSA. Annotations like sequence type and title
parsed from the PIR entry are stored in the `AnnotatedSequence`.
"""
struct PIRSequences <: SequenceFormat end

"""
    RawSequences <: SequenceFormat

A format type for reading a collection of sequences in Raw format.
Each line in the input file is treated as a separate sequence and parsed as an
individual `AnnotatedSequence`. Sequence identifiers are assigned numerically.
This contrasts with the `Raw` type (which is an `MSAFormat`), where a Raw file
is interpreted as a single multiple sequence alignment.

Use `RawSequences` when you want to process a list of sequences from a raw text file
(one sequence per line), rather than treating them as an aligned MSA.
"""
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

"""
    Utils.parse_file(io, format::Type{T}; ...) where T<:SequenceFormat -> Vector{AnnotatedSequence}

Parses a file or IO stream `io` containing a collection of sequences in a format
specified by `T` (where `T` is `FASTASequences`, `PIRSequences`, or `RawSequences`).
Returns a `Vector{AnnotatedSequence}`, where each element represents a sequence
read from the file.

This method internally uses the parsing logic of the corresponding `MSAFormat`
(e.g., `FASTA` for `FASTASequences`) to load all sequences and their potential
annotations. It then converts each sequence into an `AnnotatedSequence` object.

Keyword arguments are passed to the underlying MSA parser (`_load_sequences` and
subsequently to the format-specific parser):
- `generatemapping::Bool` (default: `false`): If true, generate "SeqMap" annotations.
- `useidcoordinates::Bool` (default: `false`): If true, parse "ID/start-end" from sequence IDs.
- `deletefullgaps::Bool` (default: `true`): If true, remove columns with 100% gaps. (Note: For `SequenceFormat`, sequences are treated individually, so this might have less impact than on an MSA).
- `keepinserts::Bool` (default: `false`): If true, try to preserve insert information (e.g., for A3M via `FASTASequences`).

```jldoctest
julia> using MIToS.MSA

julia> fasta_data = ">seq1\\nARND\\n>seq2\\nCGQH";

julia> sequences = parse_file(IOBuffer(fasta_data), FASTASequences);

julia> length(sequences)
2

julia> stringsequence(sequences[1])
"ARND"

julia> stringsequence(sequences[2])
"CGQH"
```
"""
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

"""
    Utils.print_file(io::IO, seqs::Vector{AnnotatedSequence}, format::Type{PIRSequences})

Prints a vector of `AnnotatedSequence` objects `seqs` to the stream `io` in PIR/NBRF format.
Each sequence is printed according to the PIR format rules:
1.  Header: `>{seq_type};{seq_id}`
2.  Title line.
3.  Sequence data (up to 80 chars per line), ending with `*`.

Sequence type and title are taken from the "Type" and "Title" annotations of each
`AnnotatedSequence` if available. Defaults are "XX" for type and empty for title.

```jldoctest
julia> using MIToS.MSA

julia> seq1 = AnnotatedSequence("seq1", "ARND");

julia> setannotsequence!(seq1, "Type", "P1");

julia> setannotsequence!(seq1, "Title", "Protein 1");

julia> seq2 = AnnotatedSequence("seq2", "CGQH");

julia> print_file(stdout, [seq1, seq2], PIRSequences)
>P1;seq1
Protein 1
ARND*
>XX;seq2

CGQH*
```
"""
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

"""
    Utils.print_file(io::IO, seqs::Vector{AnnotatedSequence}, format::Type{RawSequences})

Prints a vector of `AnnotatedSequence` objects `seqs` to the stream `io` in Raw format.
Each sequence is printed on a new line. No identifiers or annotations are printed.

```jldoctest
julia> using MIToS.MSA

julia> seq1 = AnnotatedSequence("seq1", "ARND");

julia> seq2 = AnnotatedSequence("seq2", "CGQH");

julia> print_file(stdout, [seq1, seq2], RawSequences)
ARND
CGQH
```
"""
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

"""
    Utils.print_file(io::IO, seqs::Vector{AnnotatedSequence}, format::Type{FASTASequences})

Prints a vector of `AnnotatedSequence` objects `seqs` to the stream `io` in FASTA format.
Each sequence is printed with a FASTA header line (`>{sequence_id}`) followed by the
sequence data on the next line(s).

```jldoctest
julia> using MIToS.MSA

julia> seq1 = AnnotatedSequence("seq1", "ARND");

julia> seq2 = AnnotatedSequence("seq2", "CGQH");

julia> print_file(stdout, [seq1, seq2], FASTASequences)
>seq1
ARND
>seq2
CGQH
```
"""
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
