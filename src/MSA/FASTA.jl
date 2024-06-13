struct FASTA <: SequenceFormat end

# FASTA Parser
# ============

function _pre_readfasta(string::AbstractString)
    seqs = split(string, '>')
    N = length(seqs) - 1
    IDS  = Array{String}(undef, N)
    SEQS = Array{String}(undef, N)
    for i in 1:N
        fields = split(seqs[i+1], '\n')
        IDS[i] = fields[1]
        SEQS[i] = replace(join(fields[2:end]), r"\s+" => "")
    end
    (IDS, SEQS)
end

# MethodError: no method matching seek(::TranscodingStreams.TranscodingStream, ::Int)
_pre_readfasta(io::TranscodingStreams.TranscodingStream) = _pre_readfasta(read(io, String))

function _pre_readfasta(io)
    IDS  = String[]
    SEQS = String[]
    for (name, seq) in FastaReader{String}(io) # FastaIO
        push!(IDS,  name)
        push!(SEQS, seq )
    end
    (IDS, SEQS)
end

function _load_sequences(io::Union{IO, AbstractString}, format::Type{FASTA}; create_annotations::Bool=false)
    IDS, SEQS = _pre_readfasta(io)
    return IDS, SEQS, Annotations()
end

# Print FASTA
# ===========

function Base.print(io::IO, msa::AbstractMatrix{Residue}, format::Type{FASTA})
    seqnames = sequencenames(msa)
    for i in 1:nsequences(msa)
        println(io, ">", seqnames[i], "\n", stringsequence(msa, i))
    end
end

Base.print(msa::AbstractMatrix{Residue}, format::Type{FASTA}) = print(stdout, msa, FASTA)
