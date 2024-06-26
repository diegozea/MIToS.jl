struct A3M <: MSAFormat end
struct A2M <: MSAFormat end

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

# A2M is similar to FASTA but uses lowercase letters and dots for inserts. In the A2M 
# format, all sequences have the same length. Since MIToS handles the inserts, we can load 
# it as FASTA. However, I will use the A3M parser instead of the FASTA parser to ensure 
# the file can be read correctly if the user confuses A2M with A3M.
_load_sequences(io::Union{IO,AbstractString}, format::Type{A2M}) = _load_sequences(io, A3M)

# Print A3M
# =========

function Utils.print_file(io::IO, msa::AbstractMatrix{Residue}, format::Union{Type{A3M},Type{A2M}})
    seqnames = sequencenames(msa)
    aligned = _get_aligned_columns(msa)
    for i in 1:nsequences(msa)
        seq = stringsequence(msa, i)
        # A2M uses dots for gaps aligned to insertions, but A3M can avoid them
        keep_insert_gaps = format === A2M
        formatted_seq = _format_inserts(seq, aligned, keep_insert_gaps)
        println(io, ">", seqnames[i])
        println(io, formatted_seq)
    end
end

Utils.print_file(msa::AbstractMatrix{Residue}, format::Union{Type{A3M},Type{A2M}}) = print_file(stdout, msa, format)
