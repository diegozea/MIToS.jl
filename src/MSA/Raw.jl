struct Raw <: MSAFormat end

# Raw Parser
# ==========

function _load_sequences(io::Union{IO, AbstractString}, format::Type{Raw}; create_annotations::Bool=false)
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

function Utils.print_file(io::IO, msa::AbstractMatrix{Residue}, format::Type{Raw})
    for i in 1:nsequences(msa)
        println(io, stringsequence(msa, i))
    end
end

Utils.print_file(msa::AbstractMatrix{Residue}, format::Type{Raw}) = print_file(stdout, msa, Raw)
