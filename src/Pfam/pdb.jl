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

"""
This function returns a `Dict{Int64,ASCIIString}` with msa column number as keys and PDB residue number as values.
The mapping is performed using SIFTS. This function needs correct *SeqMap* annotations.
If you are working with a **downloaded Pfam MSA without modifications**, you should `read` it using `generatemapping=true` and `useidcoordinates=true`.
"""
function msacolumn2pdbresidue(seqid::ASCIIString,
    pdbid::ASCIIString,
    chain::ASCIIString,
    pfamid::ASCIIString,
    msa::AnnotatedMultipleSequenceAlignment,
    siftsfile::ASCIIString)
    siftsmap = siftsmapping(siftsfile, dbPfam, pfamid, dbPDB, lowercase(pdbid), chain=chain, missings=false)
    seqmap   = getsequencemapping(msa, seqid)
    N = ncolumns(msa)
    m = Dict{Int,ASCIIString}()
    sizehint!(m, N)
    for i in 1:N
      m[i] = get(siftsmap, seqmap[i], "")
    end
    m
end

"If you don't indicate the Pfam accession number (`pfamid`), this function tries to read the *AC* file annotation."
msacolumn2pdbresidue(seqid::ASCIIString, pdbid::ASCIIString, chain::ASCIIString,
           msa::AnnotatedMultipleSequenceAlignment, siftsfile::ASCIIString) = msacolumn2pdbresidue(seqid, pdbid, chain,
           ascii(split(getannotfile(msa, "AC"), '.')[1]), msa, siftsfile::ASCIIString)

"If you don't indicate the path to the `siftsfile` used in the mapping, this function downloads the SIFTS file in the current folder."
msacolumn2pdbresidue(seqid::ASCIIString, pdbid::ASCIIString, chain::ASCIIString,
           pfamid::ASCIIString, msa::AnnotatedMultipleSequenceAlignment) = msacolumn2pdbresidue(seqid, pdbid, chain,
           pfamid, msa, downloadsifts(pdbid))

msacolumn2pdbresidue(seqid::ASCIIString, pdbid::ASCIIString, chain::ASCIIString,
           msa::AnnotatedMultipleSequenceAlignment) = msacolumn2pdbresidue(seqid, pdbid, chain,
           ascii(split(getannotfile(msa, "AC"), '.')[1]), msa)

# PDB contacts for each column
# ============================

"""
This function takes two Dicts.
The first is an `OrderedDict{ASCIIString,PDBResidue}` from PDB residue number to `PDBResidue`.
The second is a `Dict{Int,ASCIIString}` from MSA column to PDB residue number.
This returns a `PairwiseListMatrix{Float64,false}` of `0.0` and `1.0` where `1.0` indicates a residue contact
(inter residue distance less or equal to 6.05 angstroms between any heavy atom). `NaN` indicates a missing value.
"""
function msacontacts(residues::OrderedDict{ASCIIString,PDBResidue}, column2residues::Dict{Int,ASCIIString})
  contacts = PairwiseListMatrices.PairwiseListMatrix(Float64, length(column2residues), false, NaN)
  @iterateupper contacts false begin

    resnumi = get(:($column2residues), i, "")
    resnumj = get(:($column2residues), j, "")
    if resnumi != "" && resnumj != "" && haskey(:($residues), resnumi) && haskey(:($residues), resnumj)
      list[k] = Float64(:($contact)(:($residues)[resnumi], :($residues)[resnumj], 6.05))
    else
      list[k] = NaN
    end

  end
  contacts
end
