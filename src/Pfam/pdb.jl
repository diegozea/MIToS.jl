# PDB ids from Pfam sequence annotations
# ======================================

const _regex_PDB_from_GS = r"PDB;\s+(\w+)\s+(\w);\s+\w+-\w+;" # i.e.: "PDB; 2VQC A; 4-73;\n"

"Generates from a Pfam `msa` a `Dict{ASCIIString, Vector{Tuple{ASCIIString,ASCIIString}}}` with sequence ID as key and list of tuples of PDB code and chain as values."
function getseq2pdb(msa::AnnotatedMultipleSequenceAlignment)
    dict = Dict{ASCIIString, Vector{Tuple{ASCIIString,ASCIIString}}}()
    for (k, v) in getannotsequence(msa)
        id, annot = k
        # i.e.: "#=GS F112_SSV1/3-112 DR PDB; 2VQC A; 4-73;\n"
        if annot == "DR" && ismatch(_regex_PDB_from_GS, v)
            for m in eachmatch(_regex_PDB_from_GS, v)
                if haskey(dict, id)
                    push!(dict[id], (m.captures[1], m.captures[2]))
                else
                    dict[id] = Tuple{ASCIIString,ASCIIString}[ (m.captures[1], m.captures[2]) ]
                end
            end
        end
    end
    sizehint!(dict, length(dict))
end
