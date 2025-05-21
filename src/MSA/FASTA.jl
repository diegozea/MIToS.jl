"""
    FASTA <: MSAFormat

Represents the FASTA multiple sequence alignment format. In this format, each sequence
is preceded by a header line starting with `>` followed by the sequence identifier and
an optional description. The sequence itself can span multiple lines.

MIToS uses `FastaIO.jl` for efficient parsing of FASTA files. When reading, sequence
identifiers are disambiguated if necessary (e.g., if duplicates exist) by appending
suffixes like `(1)`, `(2)`, etc. The original name is stored in annotations if changed.
"""
struct FASTA <: MSAFormat end

# FASTA Parser
# ============

function _pre_readfasta(string::AbstractString)
    seqs = split(string, '>')
    N = length(seqs) - 1
    IDS = Array{String}(undef, N)
    SEQS = Array{String}(undef, N)
    for i = 1:N
        fields = split(seqs[i+1], '\n')
        IDS[i] = fields[1]
        SEQS[i] = replace(join(fields[2:end]), r"\s+" => "")
    end

    (IDS, SEQS)
end

# MethodError: no method matching seek(::TranscodingStreams.TranscodingStream, ::Int)
_pre_readfasta(io::TranscodingStreams.TranscodingStream) = _pre_readfasta(read(io, String))

function _pre_readfasta(io)
    IDS = String[]
    SEQS = String[]
    for (name, seq) in FastaReader{String}(io) # FastaIO
        push!(IDS, name)
        push!(SEQS, seq)
    end
    (IDS, SEQS)
end

function _load_sequences(
    io::Union{IO,AbstractString},
    format::Type{FASTA};
    create_annotations::Bool = false,
)
    IDS, SEQS = _pre_readfasta(io)
    annot = Annotations()
    _disambiguate_seqnames!(IDS, annot)
    return IDS, SEQS, annot
end

# Print FASTA
# ===========

"""
    Utils.print_file(io::IO, msa::AbstractMatrix{Residue}, format::Type{FASTA})

Prints the multiple sequence alignment `msa` to the stream `io` in FASTA format.
Each sequence is preceded by a header line `>` followed by its sequence identifier.
The sequence itself is printed on the next line.

```jldoctest
julia> using MIToS.MSA

julia> msa_matrix = Residue[Residue('A') Residue('R'); Residue('C') Residue('D')];

julia> named_msa = MultipleSequenceAlignment(NamedArray(msa_matrix, (["seq1", "seq2"], ["col1", "col2"])));

julia> print_file(stdout, named_msa, FASTA)
>seq1
AR
>seq2
CD
```
"""
function Utils.print_file(io::IO, msa::AbstractMatrix{Residue}, format::Type{FASTA})
    seqnames = sequencenames(msa)
    for i = 1:nsequences(msa)
        println(io, ">", seqnames[i], "\n", stringsequence(msa, i))
    end
end
