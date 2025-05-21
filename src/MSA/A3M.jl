"""
    A3M <: MSAFormat

Represents the A3M (Alignment in 3 Minutes) multiple sequence alignment format.
This format is similar to FASTA, but with specific conventions for representing
insertions relative to a consensus model, often used in HH-suite tools.

Key characteristics:
- Uppercase letters: Match states (aligned to consensus).
- Lowercase letters: Insert states (not aligned to consensus).
- `-`: Gaps in match states.
- `.`: Gaps in insert states (these are often removed or handled by `_add_insert_gaps!`).

MIToS parsing for A3M handles potential inconsistencies in sequence lengths that arise
from inserts by padding shorter sequences with gap characters (`.`) in insert regions
if necessary (`_add_insert_gaps!`). Sequence identifiers are disambiguated as in FASTA.
"""
struct A3M <: MSAFormat end

"""
    A2M <: MSAFormat

Represents the A2M multiple sequence alignment format. This format is closely related
to A3M and FASTA. Like A3M, it uses uppercase letters for match states and lowercase
letters for insert states. A key difference from A3M is that in A2M, all sequences
are expected to have the same length, with `.` characters used for gaps aligned to
insertions in other sequences.

MIToS parsing for A2M uses the same underlying mechanism as A3M, as it can correctly
handle the A2M conventions (lowercase for inserts, dots for insert gaps).
"""
struct A2M <: MSAFormat end

function _add_insert_gaps!(SEQS)
    seq_len = length.(SEQS)
    ncol = maximum(seq_len)
    j = 1
    is_insert = false
    while j <= ncol
        for i = 1:length(SEQS)
            if j <= seq_len[i]
                res = SEQS[i][j]
                if islowercase(res) || res == '.'
                    is_insert = true
                    break
                end
            end
        end
        if is_insert
            for i = 1:length(SEQS)
                seq = SEQS[i]
                if j > seq_len[i]
                    SEQS[i] = seq * "."
                    seq_len[i] += 1
                else
                    res = seq[j]
                    if isuppercase(res) || res == '-'
                        SEQS[i] = seq[1:j-1] * "." * seq[j:end]
                        seq_len[i] += 1
                    end
                end
            end
            ncol = maximum(seq_len)
        end
        is_insert = false
        j += 1
    end
    SEQS
end

function _load_sequences(
    io::Union{IO,AbstractString},
    format::Type{A3M};
    create_annotations::Bool = false,
)
    IDS, SEQS = _pre_readfasta(io)
    _check_seq_and_id_number(IDS, SEQS)
    try
        _check_seq_len(IDS, SEQS)
    catch
        SEQS = _add_insert_gaps!(SEQS)
    end
    annot = Annotations()
    _disambiguate_seqnames!(IDS, annot)
    return IDS, SEQS, annot
end

# A2M is similar to FASTA but uses lowercase letters and dots for inserts. In the A2M 
# format, all sequences have the same length. Since MIToS handles the inserts, we can load 
# it as FASTA. However, I will use the A3M parser instead of the FASTA parser to ensure 
# the file can be read correctly if the user confuses A2M with A3M.
_load_sequences(io::Union{IO,AbstractString}, format::Type{A2M}) = _load_sequences(io, A3M)

# Print A3M
# =========

"""
    Utils.print_file(io::IO, msa::AbstractMatrix{Residue}, format::Union{Type{A3M},Type{A2M}})

Prints the multiple sequence alignment `msa` to the stream `io` in A3M or A2M format.
The output is FASTA-like, with each sequence preceded by a header line `>` followed
by its sequence identifier.

The sequence formatting depends on the "Aligned" column annotation if present (typically
indicating consensus match vs. insert columns):
- Residues in "match" columns (annotated as '1' in "Aligned") are printed as uppercase.
- Residues in "insert" columns (annotated as '0' in "Aligned") are printed as lowercase.
- Gaps ('-') in "insert" columns are printed as '.' if `format` is `A2M`, otherwise they
  are omitted if `format` is `A3M` (unless `keep_insert_gaps` was true during parsing,
  which affects the `_format_inserts` behavior not explicitly shown here but implied).

```jldoctest
julia> using MIToS.MSA

# Create a sample MSA
julia> msa = AnnotatedMultipleSequenceAlignment(NamedArray(Residue["ARND"; "CGQH"]))
AnnotatedMultipleSequenceAlignment with 0 annotations : 2×4 Named Matrix{Residue}
Seq │ Col
─── │ ─────────
    │ 1   2   3   4
─── │ ─────────────
1   │ A   R   N   D
2   │ C   G   Q   H

julia> setannotcolumn!(msa, "Aligned", "1010") # A and N are match, R and D are inserts for seq1

julia> print_file(stdout, msa, A3M)
>1
ArNd
>2
CgQh

julia> print_file(stdout, msa, A2M) # Similar output, but insert gaps would be '.'
>1
ArNd
>2
CgQh
```
"""
function Utils.print_file(
    io::IO,
    msa::AbstractMatrix{Residue},
    format::Union{Type{A3M},Type{A2M}},
)
    seqnames = sequencenames(msa)
    aligned = _get_aligned_columns(msa)
    for i = 1:nsequences(msa)
        seq = stringsequence(msa, i)
        # A2M uses dots for gaps aligned to insertions, but A3M can avoid them
        keep_insert_gaps = format === A2M
        formatted_seq = _format_inserts(seq, aligned, keep_insert_gaps)
        println(io, ">", seqnames[i])
        println(io, formatted_seq)
    end
end
