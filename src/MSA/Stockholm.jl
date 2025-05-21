"""
    Stockholm <: MSAFormat

Represents the Stockholm multiple sequence alignment format. This format can store
sequence alignments along with various types of annotations for the file itself (GF),
sequences (GS), columns (GC), and residues (GR).

MIToS parses these annotations into an `Annotations` object, typically stored within
an `AnnotatedMultipleSequenceAlignment`. When printing, MIToS can also write out
these annotations in the Stockholm format. Sequences are expected to be aligned,
and may be split across multiple blocks in the file, keyed by sequence names.
MIToS handles the reassembly of such fragmented sequences.
"""
struct Stockholm <: MSAFormat end

# NOTE: Sequence‑Name Disambiguation
# We do not support sequence‑name disambiguation via the `OnlineSequenceNameDisambiguator`
# in Stockholm format, because duplicate sequence names are not permitted.
# Long sequences are split across multiple blocks, and the name is used to
# reassemble the fragments and to key per‑sequence annotations. An example of
# file with multiple blocks is at test/data/clustalo-I20240512-trunc.aln-stockholm

@inline function _fill_with_sequence_line!(IDS, SEQS, line)
    if !startswith(line, '#') && !startswith(line, "//")
        words = get_n_words(line, 2)
        @inbounds id = words[1]
        if id in IDS
            # It's useful when sequences are split into several lines
            # It can be a problem with duplicated IDs
            i = something(findfirst(isequal(id), IDS), 0)
            SEQS[i] = SEQS[i] * words[2]
        else
            push!(IDS, id)
            push!(SEQS, words[2])
        end
    end
end

function _fill_with_line!(IDS, SEQS, GF, GS, GC, GR, line)
    if startswith(line, "#=GF")
        words = get_n_words(line, 3)
        id = words[2]
        if id in keys(GF)
            GF[id] = GF[id] * "\n" * words[3]
        else
            GF[id] = words[3]
        end
    elseif startswith(line, "#=GS")
        words = get_n_words(line, 4)
        idtuple = (words[2], words[3])
        if idtuple in keys(GS)
            GS[idtuple] = GS[idtuple] * "\n" * words[4]
        else
            GS[idtuple] = words[4]
        end
    elseif startswith(line, "#=GC")
        words = get_n_words(line, 3)
        GC[words[2]] = words[3]
    elseif startswith(line, "#=GR")
        words = get_n_words(line, 4)
        GR[(words[2], words[3])] = words[4]
    else
        _fill_with_sequence_line!(IDS, SEQS, line)
    end
end

function _pre_readstockholm(io::Union{IO,AbstractString})
    IDS = OrderedSet{String}()
    SEQS = String[]
    GF = OrderedDict{String,String}()
    GC = Dict{String,String}()
    GS = Dict{Tuple{String,String},String}()
    GR = Dict{Tuple{String,String},String}()

    @inbounds for line::String in lineiterator(io)
        isempty(line) && continue
        startswith(line, "//") && break
        _fill_with_line!(IDS, SEQS, GF, GS, GC, GR, line)
    end

    GF = sizehint!(GF, length(GF))
    GC = sizehint!(GC, length(GC))
    GS = sizehint!(GS, length(GS))
    GR = sizehint!(GR, length(GR))
    (IDS, SEQS, GF, GS, GC, GR)
end

function _pre_readstockholm_sequences(io::Union{IO,AbstractString})
    IDS = OrderedSet{String}()
    SEQS = String[]
    @inbounds for line::String in lineiterator(io)
        isempty(line) && continue
        startswith(line, "//") && break
        _fill_with_sequence_line!(IDS, SEQS, line)
    end
    (IDS, SEQS)
end

function _load_sequences(
    io::Union{IO,AbstractString},
    format::Type{Stockholm};
    create_annotations::Bool = false,
)
    if create_annotations
        IDS, SEQS, GF, GS, GC, GR = _pre_readstockholm(io)
        annot = Annotations(GF, GS, GC, GR)
    else
        IDS, SEQS = _pre_readstockholm_sequences(io)
        annot = Annotations()
    end
    return collect(IDS), SEQS, annot
end

# Print Pfam
# ==========

function _to_sequence_dict(annotation::Dict{Tuple{String,String},String})
    seq_dict = Dict{String,Vector{String}}()
    for (key, value) in annotation
        seq_id = key[1]
        if haskey(seq_dict, seq_id)
            push!(seq_dict[seq_id], string(seq_id, '\t', key[2], '\t', value))
        else
            seq_dict[seq_id] = [string(seq_id, '\t', key[2], '\t', value)]
        end
    end
    sizehint!(seq_dict, length(seq_dict))
end

"""
    Utils.print_file(io::IO, msa::AbstractMatrix{Residue}, format::Type{Stockholm})

Prints the multiple sequence alignment `msa` to the stream `io` in Stockholm format.

If `msa` is an `AnnotatedAlignedObject` (e.g., `AnnotatedMultipleSequenceAlignment`),
its annotations are printed:
- File annotations (`#=GF`)
- Sequence annotations (`#=GS`)
- Residue annotations (`#=GR`), grouped by sequence.
- Column annotations (`#=GC`)

Sequences are printed with their identifiers. If an "Aligned" column annotation
exists (typically indicating insert columns from HMM profiles), residues in
insert columns are printed in lowercase, and gap characters in insert columns
can be represented as '.' if `keep_insert_gaps` was true during parsing (though this
print function currently implies `keep_insert_gaps=true` for output formatting of inserts).

The output is terminated by `//`.

```jldoctest
julia> using MIToS.MSA

julia> msa = AnnotatedMultipleSequenceAlignment(NamedArray(Residue["AR"; "CD"]))
AnnotatedMultipleSequenceAlignment with 0 annotations : 2×2 Named Matrix{Residue}
Seq │ Col
─── │ ──────
    │ 1   2
─── │ ───────
1   │ A   R
2   │ C   D

julia> setannotfile!(msa, "ID", "ExamplePfam");

julia> setannotsequence!(msa, "1", "AC", "PF00000");

julia> setannotresidue!(msa, "2", "SS", "HH");

julia> setannotcolumn!(msa, "CONS", "..");

julia> print_file(stdout, msa, Stockholm)
#=GF ID	ExamplePfam
#=GS 1	AC	PF00000
1			AR
2			CD
#=GR 2	SS	HH
#=GC CONS			..
//
```
"""
function Utils.print_file(io::IO, msa::AbstractMatrix{Residue}, format::Type{Stockholm})
    has_annotations = isa(msa, AnnotatedAlignedObject) && !isempty(msa.annotations)
    if has_annotations
        _printfileannotations(io, msa.annotations)
        _printsequencesannotations(io, msa.annotations)
        res_annotations = _to_sequence_dict(msa.annotations.residues)
    end
    seqnames = sequencenames(msa)
    aligned = _get_aligned_columns(msa)
    for i = 1:nsequences(msa)
        id = seqnames[i]
        seq = stringsequence(msa, i)
        formatted_seq = _format_inserts(seq, aligned)
        println(io, id, "\t\t\t", formatted_seq)
        if has_annotations && haskey(res_annotations, id)
            for line in res_annotations[id]
                println(io, "#=GR ", line)
            end
        end
    end
    has_annotations && _printcolumnsannotations(io, msa.annotations)
    println(io, "//")
end
