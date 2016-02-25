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
This function returns a `Dict{Int64,ASCIIString}` with **MSA column numbers on the input file** as keys and PDB residue numbers (`""` for missings) as values.
The mapping is performed using SIFTS. This function needs correct *ColMap* and *SeqMap* annotations.
This checks correspondence of the residues between the sequence and SIFTS (It throws a warning if the are differences).
If you are working with a **downloaded Pfam MSA without modifications**, you should `read` it using `generatemapping=true` and `useidcoordinates=true`.
"""
function msacolumn2pdbresidue(msa::AnnotatedMultipleSequenceAlignment,
    seqid::ASCIIString,
    pdbid::ASCIIString,
    chain::ASCIIString,
    pfamid::ASCIIString,
    siftsfile::ASCIIString)

  siftsres = read(siftsfile, SIFTSXML, chain=chain, missings=true)

  up2res = Dict{Int,Tuple{ASCIIString,ASCIIString}}()
  for res in siftsres
    if !isnull(res.Pfam) && get(res.Pfam).id == uppercase(pfamid)
      pfnum  = get(res.Pfam).number
      pfname = get(res.Pfam).name
      if !isnull(res.PDB) && (get(res.PDB).id == lowercase(pdbid)) && !res.missing
        up2res[pfnum] = (pfname, get(res.PDB).number)
      else
        up2res[pfnum] = (pfname, "")
      end
    end
  end

  seq      = ascii(getsequence(msa,  seqid))
  seqmap   = getsequencemapping(msa, seqid)
  colmap   = getcolumnmapping(msa)
  N        = ncolumns(msa)

  m = Dict{Int,ASCIIString}()
  sizehint!(m, N)
  for i in 1:N
    up_number = seqmap[i]
    if up_number != 0
      up_res, pdb_resnum = get(up2res, up_number, ("",""))
      if string(seq[i]) == up_res
        m[colmap[i]] = pdb_resnum
      else
        warn(string(pfamid, " ", seqid, " ", pdbid, " ", chain, " : MSA sequence residue at ", i, " (", seq[i], ") != SIFTS residue (UniProt/Pfam: ", up_res, ", PDB: ", pdb_resnum, ")"))
      end
    end
  end
  m
end

"If you don't indicate the path to the `siftsfile` used in the mapping, this function downloads the SIFTS file in the current folder."
msacolumn2pdbresidue(msa::AnnotatedMultipleSequenceAlignment,
                     seqid::ASCIIString, pdbid::ASCIIString, chain::ASCIIString,
                     pfamid::ASCIIString) = msacolumn2pdbresidue(msa, seqid, pdbid, chain, pfamid, downloadsifts(pdbid))

"If you don't indicate the Pfam accession number (`pfamid`), this function tries to read the *AC* file annotation."
msacolumn2pdbresidue(msa::AnnotatedMultipleSequenceAlignment,
                     seqid::ASCIIString, pdbid::ASCIIString, chain::ASCIIString) = msacolumn2pdbresidue(msa, seqid, pdbid, chain,
           ascii(split(getannotfile(msa, "AC"), '.')[1]))

"Returns a `BitVector` where there is a `true` for each column with PDB residue."
function hasresidues(msa::AnnotatedMultipleSequenceAlignment, column2residues::Dict{Int,ASCIIString})
  colmap = getcolumnmapping(msa)
  ncol = length(colmap)
  mask = falses(ncol)
  for i in 1:ncol
    if get(column2residues, colmap[i], "") != ""
      mask[i] = true
    end
  end
  mask
end

# PDB residues for each column
# ============================

"""
This function takes an `AnnotatedMultipleSequenceAlignment` with correct *ColMap* annotations and two dicts:

1. The first is an `OrderedDict{ASCIIString,PDBResidue}` from PDB residue number to `PDBResidue`.

2. The second is a `Dict{Int,ASCIIString}` from MSA column number **on the input file** to PDB residue number.

This returns an `OrderedDict{Int,PDBResidue}` from input column number (ColMap) to `PDBResidue`.
Residues on iserts are not included.
"""
function msaresidues(msa::AnnotatedMultipleSequenceAlignment, residues::OrderedDict{ASCIIString,PDBResidue}, column2residues::Dict{Int,ASCIIString})
  colmap = getcolumnmapping(msa)
  msares = sizehint!(OrderedDict{Int,PDBResidue}(), length(colmap))
  for col in colmap
    resnum = get(column2residues, col, "")
    if resnum != ""
      msares[col] = residues[resnum]
    end
  end
  sizehint!(msares, length(msares))
end

# Contact Map
# ===========

"""
This function takes an `AnnotatedMultipleSequenceAlignment` with correct *ColMap* annotations and two dicts:

1. The first is an `OrderedDict{ASCIIString,PDBResidue}` from PDB residue number to `PDBResidue`.

2. The second is a `Dict{Int,ASCIIString}` from **MSA column number on the input file** to PDB residue number.

This returns a `PairwiseListMatrix{Float64,false}` of `0.0` and `1.0` where `1.0` indicates a residue contact
(inter residue distance less or equal to 6.05 angstroms between any heavy atom). `NaN` indicates a missing value.
"""
function msacontacts(msa::AnnotatedMultipleSequenceAlignment, residues::OrderedDict{ASCIIString,PDBResidue}, column2residues::Dict{Int,ASCIIString}, distance_limit::Float64=6.05)
  colmap   = getcolumnmapping(msa)
  contacts = PairwiseListMatrices.PairwiseListMatrix(Float64, length(column2residues), colmap, false, NaN)
  @inbounds @iterateupper contacts false begin

    resnumi = get(:($column2residues), :($colmap)[i], "")
    resnumj = get(:($column2residues), :($colmap)[j], "")
    if resnumi != "" && resnumj != "" && haskey(:($residues), resnumi) && haskey(:($residues), resnumj)
      list[k] = Float64(:($contact)(:($residues)[resnumi], :($residues)[resnumj], :($distance_limit)))
    else
      list[k] = NaN
    end

  end
  contacts
end

# AUC (contact prediction)
# ========================

"""
This function takes a `msacontacts` or its list of contacts `contact_list` with 1.0 for true contacts and 0.0 for not contacts (NaN or other numbers for missing values).
Returns two `BitVector`s, the first with `true`s where `contact_list` is 1.0 and the second with `true`s where `contact_list` is 0.0. There are useful for AUC calculations.
"""
function getcontactmasks{T <: AbstractFloat}(contact_list::Vector{T})
  N = length(contact_list)
  true_contacts  = falses(N)
  false_contacts = trues(N) # In general, there are few contacts
  @inbounds for i in 1:N
    value = contact_list[i]
    if value == 1.0
      true_contacts[i] = true
    end
    if value != 0.0
      false_contacts[i] = false
    end
  end
  true_contacts, false_contacts
end

getcontactmasks{T <: AbstractFloat}(msacontacts::PairwiseListMatrix{T,false}) = getcontactmasks(getlist(msacontacts))

"""
Returns the Area Under a ROC (Receiver Operating Characteristic) Curve (AUC) of the `scores_list` for `true_contacts` prediction.
The three vectors should have the same length and `false_contacts` should be `true` where there are not contacts.
"""
AUC{T}(scores_list::Vector{T}, true_contacts::BitVector, false_contacts::BitVector) = 1 - auc(roc(scores_list[true_contacts  & !isnan(scores_list)],
                                                                                                  scores_list[false_contacts & !isnan(scores_list)]))

"""
Returns the Area Under a ROC (Receiver Operating Characteristic) Curve (AUC) of the `scores` for `true_contacts` prediction.
`scores`, `true_contacts` and `false_contacts` should have the same number of elements and `false_contacts` should be `true` where there are not contacts.
"""
AUC{T}(scores::PairwiseListMatrix{T,false}, true_contacts::BitVector, false_contacts::BitVector) = AUC(getlist(scores), true_contacts, false_contacts)

"""
Returns the Area Under a ROC (Receiver Operating Characteristic) Curve (AUC) of the `scores` for `msacontact` prediction.
`score` and `msacontact` lists are vinculated (inner join) by their labels (i.e. column number in the file).
`msacontact` should have 1.0 for true contacts and 0.0 for not contacts (NaN or other numbers for missing values).
"""
function AUC{S <: AbstractFloat, C <: AbstractFloat}(scores::PairwiseListMatrix{S,false}, msacontacts::PairwiseListMatrix{C,false})
  sco, con = join(scores, msacontacts, kind=:inner)
  true_contacts, false_contacts = getcontactmasks(con)
  AUC(sco, true_contacts, false_contacts)
end
