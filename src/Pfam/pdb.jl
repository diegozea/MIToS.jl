# PDB ids from Pfam sequence annotations
# ======================================

const _regex_PDB_from_GS = r"PDB;\s+(\w+)\s+(\w);\s+\w+-\w+;" # i.e.: "PDB; 2VQC A; 4-73;\n"

"""Generates from a Pfam `msa` a `Dict{ASCIIString, Vector{Tuple{ASCIIString,ASCIIString}}}`.
Keys are sequence IDs and each value is a list of tuples containing PDB code and chain.

```
julia> getseq2pdb(msa)
Dict{ASCIIString,Array{Tuple{ASCIIString,ASCIIString},1}} with 1 entry:
  "F112_SSV1/3-112" => [("2VQC","A")]

```
"""
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

# Mapping PDB/Pfam
# ================

function getcol2res(seqid::ASCIIString,
    pdbid::ASCIIString,
    chain::ASCIIString,
    pfamid::ASCIIString,
    msa::AnnotatedMultipleSequenceAlignment,
    siftsfile::ASCIIString)
    siftsmap = siftsmapping(siftsfile, dbPfam, pfamid, dbPDB, lowercase(pdbid), chain=chain, missings=false)
    seqmap   = getsequencemapping(msa, seqid)
    colmap   = getcolumnmapping(msa)
    N = length(colmap)
    m = Dict{Int,ASCIIString}()
    sizehint!(m, N)
    for i in 1:N
      m[colmap[i]] = get(siftsmap, seqmap[i], "")
    end
    m
end

getcol2res(seqid::ASCIIString, pdbid::ASCIIString, chain::ASCIIString,
           msa::AnnotatedMultipleSequenceAlignment, siftsfile::ASCIIString) = getcol2res(seqid, pdbid, chain,
           ascii(split(getannotfile(msa, "AC"), '.')[1]), msa, siftsfile::ASCIIString)

getcol2res(seqid::ASCIIString, pdbid::ASCIIString, chain::ASCIIString,
           pfamid::ASCIIString, msa::AnnotatedMultipleSequenceAlignment) = getcol2res(seqid, pdbid, chain,
           pfamid, msa, downloadsifts(pdbid))

getcol2res(seqid::ASCIIString, pdbid::ASCIIString, chain::ASCIIString,
           msa::AnnotatedMultipleSequenceAlignment) = getcol2res(seqid, pdbid, chain,
           ascii(split(getannotfile(msa, "AC"), '.')[1]), msa)

