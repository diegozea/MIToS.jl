struct FASTA <: Format end

# FASTA Parser
# ============

function _pre_readfasta(io::AbstractString)
    seqs = split(io, '>')
    N = length(seqs) - 1
    IDS  = Array{String}(N)
    SEQS = Array{String}(N)
    for i in 1:N
        fields = split(seqs[i+1], '\n')
        IDS[i] = fields[1]
        SEQS[i] = replace(fields[2], r"\s+", "")
    end
    _check_seq_len(IDS, SEQS)
    (IDS, SEQS)
end

function _pre_readfasta(io::IO)
    IDS  = String[]
    SEQS = String[]
    for (name, seq) in FastaReader{String}(io) # FastaIO
        push!(IDS,  name)
        push!(SEQS, seq )
    end
    _check_seq_len(IDS, SEQS)
    (IDS, SEQS)
end

function Base.parse(io::Union{IO, AbstractString},
                    format::Type{FASTA},
                    output::Type{AnnotatedMultipleSequenceAlignment};
                    generatemapping::Bool=false,
                    useidcoordinates::Bool=false,
                    deletefullgaps::Bool=true,
                    keepinserts::Bool=false)
    IDS, SEQS = _pre_readfasta(io)
    annot = Annotations()
    _generate_annotated_msa(annot, IDS, SEQS, keepinserts, generatemapping,
        useidcoordinates, deletefullgaps)
end

function Base.parse(io::Union{IO, AbstractString},
                    format::Type{FASTA},
                    output::Type{NamedResidueMatrix};
                    deletefullgaps::Bool=true)
    IDS, SEQS = _pre_readfasta(io)
    msa = _generate_named_array(SEQS, IDS)
    if deletefullgaps
        return deletefullgapcolumns(msa)
    end
    msa
end

function Base.parse(io::Union{IO, AbstractString},
                   format::Type{FASTA},
                   output::Type{MultipleSequenceAlignment};
                   deletefullgaps::Bool=true)
    msa = parse(io, format, NamedResidueMatrix, deletefullgaps=deletefullgaps)
    MultipleSequenceAlignment(msa)
end

function Base.parse(io::Union{IO, AbstractString},
                    format::Type{FASTA},
                    output::Type{Matrix{Residue}};
                    deletefullgaps::Bool=true)
    IDS, SEQS = _pre_readfasta(io)
    _strings_to_matrix_residue_unsafe(SEQS, deletefullgaps)
end

function Base.parse(io::Union{IO, AbstractString},
                    format::Type{FASTA};
                    generatemapping::Bool=false,
                    useidcoordinates::Bool=false,
                    deletefullgaps::Bool=true,
                    keepinserts::Bool=false)
    parse(io, FASTA, AnnotatedMultipleSequenceAlignment; generatemapping=generatemapping,
        useidcoordinates=useidcoordinates, deletefullgaps=deletefullgaps,
        keepinserts=keepinserts)
end

# Print FASTA
# ===========

function Base.print(io::IO, msa::AbstractMatrix{Residue}, format::Type{FASTA})
    seqnames = sequencenames(msa)
    for i in 1:nsequences(msa)
        println(io, string(">", seqnames[i], "\n", stringsequence(msa, i)))
    end
end

Base.print(msa::AbstractMatrix{Residue}, format::Type{FASTA}) = print(STDOUT, msa, FASTA)
