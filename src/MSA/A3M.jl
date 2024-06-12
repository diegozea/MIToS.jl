struct A3M <: FileFormat end

function _add_inserts(SEQS)
    seq_len = length.(SEQS)
    ncol = maximum(seq_len)
    j = 1
    is_insert = false
    while j <= ncol
        for i in 1:length(SEQS)
            if j <= seq_len[i]
                res = SEQS[i][j]
                if islowercase(res) || res == '.'
                    is_insert = true
                    break
                end
            end
        end
        if is_insert
            for i in 1:length(SEQS)
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

function _pre_reada3m(io)
    IDS, SEQS = _pre_readfasta(io)
    _check_seq_and_id_number(IDS, SEQS)
    try
        _check_seq_len(IDS, SEQS)
    catch
        SEQS = _add_inserts(SEQS)
    end
    return IDS, SEQS
end

# ----------------- the following functions are copied from FASTA.jl -----------------

# TODO: refactor to use multiple dispatch to avoid code duplication

function Base.parse(io::Union{IO, AbstractString},
    format::Type{A3M},
    output::Type{AnnotatedMultipleSequenceAlignment};
    generatemapping::Bool=false,
    useidcoordinates::Bool=false,
    deletefullgaps::Bool=true,
    keepinserts::Bool=false)
IDS, SEQS = _pre_reada3m(io)
_check_seq_len(IDS, SEQS)
annot = Annotations()
_generate_annotated_msa(annot, IDS, SEQS, keepinserts, generatemapping,
useidcoordinates, deletefullgaps)
end

function Base.parse(io::Union{IO, AbstractString},
    format::Type{A3M},
    output::Type{NamedResidueMatrix{Array{Residue,2}}};
    deletefullgaps::Bool=true)
IDS, SEQS = _pre_reada3m(io)
_check_seq_len(IDS, SEQS)
msa = _generate_named_array(SEQS, IDS)
if deletefullgaps
return deletefullgapcolumns(msa)
end
msa
end

function Base.parse(io::Union{IO, AbstractString},
   format::Type{A3M},
   output::Type{MultipleSequenceAlignment};
   deletefullgaps::Bool=true)
msa = parse(io, format, NamedResidueMatrix{Array{Residue,2}},
deletefullgaps=deletefullgaps)
MultipleSequenceAlignment(msa)
end

function Base.parse(io::Union{IO, AbstractString},
    format::Type{A3M},
    output::Type{Matrix{Residue}};
    deletefullgaps::Bool=true)
IDS, SEQS = _pre_reada3m(io)
_check_seq_len(IDS, SEQS)
_strings_to_matrix_residue_unsafe(SEQS, deletefullgaps)
end

function Base.parse(io::Union{IO, AbstractString},
    format::Type{A3M};
    generatemapping::Bool=false,
    useidcoordinates::Bool=false,
    deletefullgaps::Bool=true,
    keepinserts::Bool=false)
parse(io, A3M, AnnotatedMultipleSequenceAlignment; generatemapping=generatemapping,
useidcoordinates=useidcoordinates, deletefullgaps=deletefullgaps,
keepinserts=keepinserts)
end

# Print A3M
# =========

function Base.print(io::IO, msa::AbstractMatrix{Residue}, format::Type{A3M})
seqnames = sequencenames(msa)
for i in 1:nsequences(msa)
println(io, ">", seqnames[i], "\n", stringsequence(msa, i))
end
end

Base.print(msa::AbstractMatrix{Residue}, format::Type{A3M}) = print(stdout, msa, A3M)
