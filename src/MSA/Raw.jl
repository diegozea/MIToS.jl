"""
    Raw <: MSAFormat

Represents a raw multiple sequence alignment format. In this format, each line of the
file is treated as a single sequence. No metadata or sequence identifiers are expected
within the file itself; sequence identifiers are typically assigned numerically during parsing.

This format is useful for simple, unannotated alignments where sequences are provided
one per line.
"""
struct Raw <: MSAFormat end

# Raw Parser
# ==========

function _load_sequences(
    io::Union{IO,AbstractString},
    format::Type{Raw};
    create_annotations::Bool = false,
)
    SEQS = String[]
    IDS = String[]
    for (i, line::String) in enumerate(lineiterator(io))
        push!(SEQS, line)
        push!(IDS, string(i))
    end
    return IDS, SEQS, Annotations()
end

# Print Raw
# =========

"""
    Utils.print_file(io::IO, msa::AbstractMatrix{Residue}, format::Type{Raw})

Prints the multiple sequence alignment `msa` to the stream `io` in Raw format.
Each sequence in the MSA is printed on a new line. No identifiers or annotations
are printed.

```jldoctest
julia> using MIToS.MSA

julia> msa_matrix = Residue[Residue('A') Residue('R'); Residue('C') Residue('D')];

julia> print_file(stdout, msa_matrix, Raw)
AR
CD
```
"""
function Utils.print_file(io::IO, msa::AbstractMatrix{Residue}, format::Type{Raw})
    for i = 1:nsequences(msa)
        println(io, stringsequence(msa, i))
    end
end
