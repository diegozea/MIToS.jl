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
    for (i, original_line::String) in enumerate(lineiterator(io))
        processed_line = strip(original_line) # Step 1: Remove leading/trailing whitespace
        if startswith(processed_line, "#")
            processed_line = processed_line[2:end] # Step 2: Remove '#' if it's at the start
        end
        processed_line = replace(processed_line, " " => "") # Step 3: Remove all spaces
        push!(SEQS, processed_line) # Step 4: Push the processed line
        push!(IDS, string(i))
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
