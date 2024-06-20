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

function Base.print(io::IO, msa::AbstractMatrix{Residue}, format::Type{Raw})
    for i in 1:nsequences(msa)
        println(io, stringsequence(msa, i))
    end
end

Base.print(msa::AbstractMatrix{Residue}, format::Type{Raw}) = print(stdout, msa, Raw)
