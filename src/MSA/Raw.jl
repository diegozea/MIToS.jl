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
    nonres = r"[^A-Za-z-\.]+"
    for (i, line) in enumerate(lineiterator(io))
        clean = replace(line, nonres => "")
        if !isempty(clean)
            push!(SEQS, clean)
            push!(IDS, string(i))
        end
    end
    return IDS, SEQS, Annotations()
end

# Print Raw
# =========

function Utils.print_file(io::IO, msa::AbstractMatrix{Residue}, format::Type{Raw})
    for i = 1:nsequences(msa)
        println(io, stringsequence(msa, i))
    end
end
