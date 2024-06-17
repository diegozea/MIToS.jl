struct A3M <: SequenceFormat end

function _add_insert_gaps!(SEQS)
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

function _load_sequences(io::Union{IO,AbstractString}, format::Type{A3M}; create_annotations::Bool=false)
    IDS, SEQS = _pre_readfasta(io)
    _check_seq_and_id_number(IDS, SEQS)
    try
        _check_seq_len(IDS, SEQS)
    catch
        SEQS = _add_insert_gaps!(SEQS)
    end
    return IDS, SEQS, Annotations()
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
