# PDB ids from Pfam sequence annotations
# ======================================

const _regex_PDB_from_GS = r"PDB;\s+(\w+)\s+(\w);\s+\w+-\w+;" # i.e.: "PDB; 2VQC A; 4-73;\n"

"""
Generates from a Pfam `msa` a `Dict{String, Vector{Tuple{String,String}}}`.
Keys are sequence IDs and each value is a list of tuples containing PDB code and chain.

```julia
julia> getseq2pdb(msa)
Dict{String,Array{Tuple{String,String},1}} with 1 entry:
  "F112_SSV1/3-112" => [("2VQC","A")]

```
"""
function getseq2pdb(msa::AnnotatedMultipleSequenceAlignment)
    dict = Dict{String,Vector{Tuple{String,String}}}()
    for (k, v) in getannotsequence(msa)
        id, annot = k
        # i.e.: "#=GS F112_SSV1/3-112 DR PDB; 2VQC A; 4-73;\n"
        if annot == "DR" && occursin(_regex_PDB_from_GS, v)
            for m in eachmatch(_regex_PDB_from_GS, v)
                if haskey(dict, id)
                    push!(dict[id], (m.captures[1], m.captures[2]))
                else
                    dict[id] = Tuple{String,String}[ (m.captures[1], m.captures[2]) ]
                end
            end
        end
    end
    sizehint!(dict, length(dict))
end

# Mapping PDB/Pfam
# ================

"""
`msacolumn2pdbresidue(msa, seqid, pdbid, chain, pfamid, siftsfile; strict=false, checkpdbname=false, missings=true)`

This function returns a `Dict{Int,String}` with **MSA column numbers on the input file**
as keys and PDB residue numbers (`""` for missings) as values. The mapping is performed
using SIFTS. This function needs correct *ColMap* and *SeqMap* annotations. This checks
correspondence of the residues between the MSA sequence and SIFTS
(It throws a warning if there are differences). Missing residues are included if the
keyword argument `missings` is `true` (default: `true`). If the keyword argument `strict`
is `true` (default: `false`), throws an Error, instead of a Warning, when residues don't
match. If the keyword argument `checkpdbname` is `true` (default: `false`), throws an Error
if the three letter name of the PDB residue isn't the MSA residue. If you are working with
a **downloaded Pfam MSA without modifications**, you should `read` it using
`generatemapping=true` and `useidcoordinates=true`. If you don't indicate the path to the
`siftsfile` used in the mapping, this function downloads the SIFTS file in the current
folder. If you don't indicate the Pfam accession number (`pfamid`), this function tries to
read the *AC* file annotation.
"""
function msacolumn2pdbresidue(msa::AnnotatedMultipleSequenceAlignment,
                              seqid::String,
                              pdbid::String,
                              chain::String,
                              pfamid::String,
                              siftsfile::String;
                              strict::Bool=false,
                              checkpdbname::Bool=false,
                              missings::Bool=true)

    siftsres = read(siftsfile, SIFTSXML, chain=chain, missings=missings)

    up2res = Dict{String,Tuple{String,String,Char}}()
    for res in siftsres
        if !ismissing(res.Pfam) && res.Pfam.id == uppercase(pfamid)
            pfnum  = res.Pfam.number
            if pfnum == ""
                continue
            end
            pfname = res.Pfam.name
            if !ismissing(res.PDB) && (res.PDB.id == lowercase(pdbid)) && !res.missing
                up2res[pfnum] = checkpdbname ?
                    (pfname,res.PDB.number,three2residue(res.PDB.name)) :
                    (pfname,res.PDB.number,'-')
            else
                up2res[pfnum] = checkpdbname ?
                    (pfname,"",ismissing(res.PDB) ? "" : three2residue(res.PDB.name)) :
                    (pfname,"",'-')
            end
        end
    end

    seq      = Char[x for x in vec(getsequence(msa, seqid))]
    seqmap   = getsequencemapping(msa, seqid)
    colmap   = getcolumnmapping(msa)
    N        = ncolumns(msa)

    m = Dict{Int,String}()
    sizehint!(m, N)
    for i in 1:N
        up_number = string(seqmap[i])
        if up_number != "0"
            up_res, pdb_resnum, pdb_res = get(up2res, up_number, ("","",'-'))
            if string(seq[i]) == up_res
                m[colmap[i]] = pdb_resnum
            else
                msg = string(pfamid, " ", seqid, " ", pdbid, " ", chain,
                             " : MSA sequence residue at ", i, " (", seq[i],
                             ") != SIFTS residue (UniProt/Pfam: ", up_res, ", PDB: ",
                             pdb_resnum, ")")
                strict ? throw(ErrorException(msg)) : warn(msg)
            end
            if ( checkpdbname && (seq[i] != pdb_res) )
                msg = string(pfamid, " ", seqid, " ", pdbid, " ", chain,
                             " : MSA sequence residue at ", i, " (", seq[i],
                             ") != PDB residue at ", pdb_resnum, " (", pdb_res, ")")
                throw(ErrorException(msg))
            end
        end
    end
    m
end

function msacolumn2pdbresidue(msa::AnnotatedMultipleSequenceAlignment,
                              seqid::String, pdbid::String, chain::String,
                              pfamid::String; kargs...)
    msacolumn2pdbresidue(msa, seqid, pdbid, chain, pfamid, downloadsifts(pdbid), kargs...)
end

function msacolumn2pdbresidue(msa::AnnotatedMultipleSequenceAlignment,
                              seqid::String, pdbid::String, chain::String; kargs...)
    msacolumn2pdbresidue(msa,seqid,pdbid,chain,
                         String(split(getannotfile(msa,"AC"),'.')[1]),
                         kargs...)
end

"Returns a `BitVector` where there is a `true` for each column with PDB residue."
function hasresidues(msa::AnnotatedMultipleSequenceAlignment,
                    column2residues::Dict{Int,String})
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
This function takes an `AnnotatedMultipleSequenceAlignment` with correct *ColMap*
annotations and two dicts:

1. The first is an `OrderedDict{String,PDBResidue}` from PDB residue number to `PDBResidue`.
2. The second is a `Dict{Int,String}` from MSA column number **on the input file** to PDB residue number.

`msaresidues` returns an `OrderedDict{Int,PDBResidue}` from input column number (ColMap)
to `PDBResidue`. Residues on inserts are not included.
"""
function msaresidues(msa::AnnotatedMultipleSequenceAlignment,
                     residues::OrderedDict{String,PDBResidue},
                     column2residues::Dict{Int,String})
    colmap = getcolumnmapping(msa)
    msares = sizehint!(OrderedDict{Int,PDBResidue}(), length(colmap))
    for col in colmap
        resnum = get(column2residues, col, "")
        if resnum != ""
            if haskey(residues, resnum)
                msares[col] = residues[resnum]
            else
                @warn("MSA column $col : The residue number $resnum isn't in the residues Dict.")
            end
        end
    end
    sizehint!(msares, length(msares))
end

# Contact Map
# ===========

"""
This function takes an `AnnotatedMultipleSequenceAlignment` with correct *ColMap*
annotations and two dicts:

1. The first is an `OrderedDict{String,PDBResidue}` from PDB residue number to `PDBResidue`.
2. The second is a `Dict{Int,String}` from **MSA column number on the input file** to PDB residue number.

`msacontacts` returns a `PairwiseListMatrix{Float64,false}` of `0.0` and `1.0` where `1.0`
indicates a residue contact. Contacts are defined with an inter residue distance less or
equal to `distance_limit` (default to `6.05`) angstroms between any heavy atom. `NaN`
indicates a missing value.
"""
function msacontacts(msa::AnnotatedMultipleSequenceAlignment,
                     residues::OrderedDict{String,PDBResidue},
                     column2residues::Dict{Int,String},
                     distance_limit::Float64 = 6.05)
    colmap   = getcolumnmapping(msa)
    contacts = columnpairsmatrix(msa)
    plm = getarray(contacts)
    @inbounds @iterateupper plm false begin

        resi = get(:($column2residues), :($colmap)[i], "")
        resj = get(:($column2residues), :($colmap)[j], "")
        if resi != "" && resj != "" && haskey(:($residues),resi) && haskey(:($residues),resj)
            list[k] = Float64(:($contact)(:($residues)[resi],:($residues)[resj],:($distance_limit)))
        else
            list[k] = NaN
        end

    end
    contacts
end

# AUC (contact prediction)
# ========================

"""
This function takes a `msacontacts` or its list of contacts `contact_list` with 1.0 for
true contacts and 0.0 for not contacts (NaN or other numbers for missing values).
Returns two `BitVector`s, the first with `true`s where `contact_list` is 1.0 and the second
with `true`s where `contact_list` is 0.0. There are useful for AUC calculations.
"""
function getcontactmasks(contact_list::Vector{T}) where T <: AbstractFloat
    N = length(contact_list)
    true_contacts  = falses(N)
    false_contacts = falses(N)
    @inbounds for i in 1:N
        value = contact_list[i]
        if value == 1.0
            true_contacts[i]  = true
        elseif value == 0.0
            false_contacts[i] = true
        end
        # If value is NaN, It keeps the false value
    end
    true_contacts, false_contacts
end

function getcontactmasks(plm::PairwiseListMatrix{T,false,VT}) where {T <: AbstractFloat,VT}
    getcontactmasks(getlist(plm))
end

function getcontactmasks(nplm::NamedArray{T,2,PairwiseListMatrix{T,false,TV},DN}) where {T,TV,DN}
    getcontactmasks(getarray(nplm))
end
