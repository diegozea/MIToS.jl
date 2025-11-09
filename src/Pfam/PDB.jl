# PDB ids from Pfam sequence annotations
# ======================================

const _regex_PDB_from_GS = r"PDB;\s+(\w+)\s+(\w);\s+\w+-\w+;" # i.e.: "PDB; 2VQC A; 4-73;\n"

"""
    getseq2pdb(
        msa::AnnotatedMultipleSequenceAlignment;
        force_sifts_mapping::Bool = false,
        sifts_pfam_csv::Union{Nothing,AbstractString} = nothing,
        sifts_uniprot_csv::Union{Nothing,AbstractString} = nothing,
    )

Build a `Dict{String, Vector{Tuple{String,String}}}` that maps each Pfam MSA sequence
name to the `(pdb_id, chain_id)` pairs annotated for that sequence.

The function first inspects the Pfam `#=GS` annotations (`DR PDB; …`) and records any
PDB/chain pairs it finds. If the alignment lacks those annotations (common in recent
Pfam releases) or if `force_sifts_mapping` is `true`, the mapping is complemented with
data from PDBe SIFTS summary CSV files:

  - `pdb_chain_pfam.csv.gz`
  - `pdb_chain_uniprot.csv.gz`

When `sifts_pfam_csv`/`sifts_uniprot_csv` are not provided, the corresponding files are
downloaded (or reused if already present in the working directory). The MSA must contain
an `AC` record in the file annotations so the Pfam accession number can be inferred;
otherwise an error is thrown. SIFTS-derived mappings are only added when the UniProt
residue ranges overlap with the ranges encoded in the Pfam sequence
names (e.g. `"ICSA_SHIFL/612-720"`).

SIFTS summary CSVs are parsed at most once per file absolute path and modification time
pair during the lifetime of the Julia session, so repeated `getseq2pdb` calls reuse
in-memory tables unless the underlying files change.

```julia
julia> getseq2pdb(msa) # PF18883
Dict{String, Vector{Tuple{String, String}}} with 1 entry:
  "ICSA_SHIFL/612-720" => [("3ML3", "A"), ("5KE1", "A"), ("5KE1", "B")]
```
"""
function getseq2pdb(
    msa::AnnotatedMultipleSequenceAlignment;
    force_sifts_mapping::Bool = false,
    sifts_pfam_csv::Union{Nothing,AbstractString} = nothing,
    sifts_uniprot_csv::Union{Nothing,AbstractString} = nothing,
)
    # Build the dict from sequence annotations first
    has_annot = false
    dict = Dict{String,Vector{Tuple{String,String}}}()
    for (k, v) in getannotsequence(msa)
        id, annot = k
        # i.e.: "#=GS F112_SSV1/3-112 DR PDB; 2VQC A; 4-73;\n"
        if annot == "DR" && occursin(_regex_PDB_from_GS, v)
            has_annot = true
            for m in eachmatch(_regex_PDB_from_GS, v)
                if haskey(dict, id)
                    push!(dict[id], (m.captures[1], m.captures[2]))
                else
                    dict[id] = Tuple{String,String}[(m.captures[1], m.captures[2])]
                end
            end
        end
    end
    if !has_annot
        first_seq = first(sequencename_iterator(msa))
        @warn "The MSA lacks DR PDB annotations (common in recent Pfam)." maxlog = 1 _id = "pdb_annot_$(first_seq)"
    end
    # If no annotations or forced, use SIFTS mapping
    if force_sifts_mapping || !has_annot
        _sifts_seq2pdb!(dict, msa, sifts_pfam_csv, sifts_uniprot_csv)
    end
    sizehint!(dict, length(dict))
end

# Same helper functions as in src/SIFTS/Summary.jl
_summary_name(db::Type{dbPfam}) = "pdb_chain_pfam.csv.gz"
_summary_name(db::Type{dbUniProt}) = "pdb_chain_uniprot.csv.gz"

function _download_or_reuse_sifts_file(
    sifts_csv::Union{Nothing,AbstractString},
    db::Type{T},
) where {T<:DataBase}
    if !isnothing(sifts_csv)
        return sifts_csv
    end
    filename = _summary_name(db)
    if isfile(filename) && filesize(filename) > 0
        @info "Reusing existing SIFTS file: $filename" maxlog = 1 _id = "reuse_$filename"
        filename
    else
        @info "Downloading SIFTS file: $filename" maxlog = 1 _id = "download_$filename"
        downloadsifts(db, filename = filename)
    end
end

# Cache SIFTS CSV tables so we only parse them once per file/mtime pair
const _SIFTS_TABLE_CACHE = Dict{
    Tuple{String,Float64},
    NamedTuple{(:colnames, :table),Tuple{Vector{Symbol},Matrix{String}}},
}()

function _read_sifts_table(path::AbstractString)
    abs_path = abspath(path)
    mod_time = mtime(abs_path)
    key = (abs_path, mod_time)
    if haskey(_SIFTS_TABLE_CACHE, key)
        return _SIFTS_TABLE_CACHE[key]
    end
    table = read_file(abs_path, SIFTSCSV)
    # Remove old entries for the same location
    filter!(kv -> kv.first[1] != abs_path, _SIFTS_TABLE_CACHE)
    _SIFTS_TABLE_CACHE[key] = table
    table
end

# Since the SIFTS mapping uses primary (citable) accession numbers from UniProt 
# (without version numbers), but the Pfam MSA uses the entry names, we need to get the
# mapping between these two identifiers. Pfam MSAs include the UniProt accession (with
# version numbers) in the sequence annotations, so we can extract it from there:
function _get_acc2seqnames(msa::AnnotatedMultipleSequenceAlignment)
    seq_annot = getannotsequence(msa)
    acc2seqnames = Dict{String,Vector{String}}()
    for ((seqname, annot), value) in seq_annot
        if annot == "AC"
            uniprot_acc = first(split(value, '.')) # Remove version number
            push!(get!(acc2seqnames, uniprot_acc, String[]), seqname)
        end
    end
    sizehint!(acc2seqnames, length(acc2seqnames))
    acc2seqnames
end

function _sifts_seq2pdb!(
    dict::Dict{String,Vector{Tuple{String,String}}},
    msa::AnnotatedMultipleSequenceAlignment,
    sifts_pfam_csv::Union{Nothing,AbstractString},
    sifts_uniprot_csv::Union{Nothing,AbstractString},
)
    # Get the Pfam accession from the MSA annotations
    msa_ac = getannotfile(msa, "AC", "")
    if isempty(msa_ac)
        throw(
            ErrorException("Cannot determine Pfam accession number from MSA annotations."),
        )
    end
    pfam_id = String(first(split(msa_ac, '.')))
    # Download or reuse SIFTS files
    sifts_pfam_file = _download_or_reuse_sifts_file(sifts_pfam_csv, dbPfam)
    sifts_uniprot_file = _download_or_reuse_sifts_file(sifts_uniprot_csv, dbUniProt)
    # Read SIFTS files
    sifts_pfam = _read_sifts_table(sifts_pfam_file)
    sifts_uniprot = _read_sifts_table(sifts_uniprot_file)
    # Get the PDB and UniProt accessions associated to the given Pfam
    pdb_chain_up_pfam = sifts_pfam.table[sifts_pfam.table[:, 4].==pfam_id, 1:4]
    # Keep only UniProt–PDB chain mappings that belong to the given Pfam
    uniprot_accs = Set(pdb_chain_up_pfam[:, 3])
    row_selector = in.(sifts_uniprot.table[:, 3], Ref{Set{String}}(uniprot_accs))
    pdb_chain_up_coords = sifts_uniprot.table[row_selector, :]
    # Create a mapping from UniProt primary accession number to sequence names
    acc2seqnames = _get_acc2seqnames(msa)
    # Fill in the dict with PDB–chain tuples from SIFTS
    for row_index in axes(pdb_chain_up_coords, 1)
        uniprot_acc = pdb_chain_up_coords[row_index, 3]
        if haskey(acc2seqnames, uniprot_acc)
            pdb_id = uppercase(pdb_chain_up_coords[row_index, 1])
            chain_id = uppercase(pdb_chain_up_coords[row_index, 2])
            for seqname in acc2seqnames[uniprot_acc]
                # Extract UniProt region of the given MSA sequence
                start_str, end_str = split(last(split(seqname, '/')), '-')
                seq_start = parse(Int, start_str)
                seq_end = parse(Int, end_str)
                # Extract UniProt region of the SIFTS mapping
                up_start = parse(Int, pdb_chain_up_coords[row_index, 8])
                up_end = parse(Int, pdb_chain_up_coords[row_index, 9])
                # Only add the mapping if there is overlap between both regions
                if (seq_end < up_start) || (seq_start > up_end)
                    continue # no overlap
                end
                if haskey(dict, seqname)
                    push!(dict[seqname], (pdb_id, chain_id))
                else
                    dict[seqname] = Tuple{String,String}[(pdb_id, chain_id)]
                end
            end
        end
    end
    dict
end

# Mapping PDB/Pfam
# ================

"""
`msacolumn2pdbresidue(msa, seqid, pdbid, chain, pfamid, siftsfile; strict=false, checkpdbname=false, missings=true)`

This function returns a `OrderedDict{Int,String}` with **MSA column numbers on the input file**
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
function msacolumn2pdbresidue(
    msa::AnnotatedMultipleSequenceAlignment,
    seqid::String,
    pdbid::String,
    chain::String,
    pfamid::String,
    siftsfile::String;
    strict::Bool = false,
    checkpdbname::Bool = false,
    missings::Bool = true,
)

    siftsres = read_file(siftsfile, SIFTSXML, chain = chain, missings = missings)

    up2res = OrderedDict{String,Tuple{String,String,Char}}()
    for res in siftsres
        if !ismissing(res.Pfam) && res.Pfam.id == uppercase(pfamid)
            pfnum = res.Pfam.number
            if pfnum == ""
                continue
            end
            pfname = res.Pfam.name
            if !ismissing(res.PDB) && (res.PDB.id == lowercase(pdbid)) && !res.missing
                up2res[pfnum] =
                    checkpdbname ? (pfname, res.PDB.number, three2residue(res.PDB.name)) :
                    (pfname, res.PDB.number, '-')
            else
                up2res[pfnum] =
                    checkpdbname ?
                    (pfname, "", ismissing(res.PDB) ? "" : three2residue(res.PDB.name)) :
                    (pfname, "", '-')
            end
        end
    end

    seq = Char[x for x in vec(getsequence(msa, seqid))]
    seqmap = getsequencemapping(msa, seqid)
    colmap = getcolumnmapping(msa)
    N = ncolumns(msa)

    m = OrderedDict{Int,String}()
    sizehint!(m, N)
    for i = 1:N
        up_number = string(seqmap[i])
        if up_number != "0"
            up_res, pdb_resnum, pdb_res = get(up2res, up_number, ("", "", '-'))
            if string(seq[i]) == up_res
                m[colmap[i]] = pdb_resnum
            else
                msg = string(
                    pfamid,
                    " ",
                    seqid,
                    " ",
                    pdbid,
                    " ",
                    chain,
                    " : MSA sequence residue at ",
                    i,
                    " (",
                    seq[i],
                    ") != SIFTS residue (UniProt/Pfam: ",
                    up_res,
                    ", PDB: ",
                    pdb_resnum,
                    ")",
                )
                strict ? throw(ErrorException(msg)) : @warn(msg)
            end
            if (checkpdbname && (seq[i] != pdb_res))
                msg = string(
                    pfamid,
                    " ",
                    seqid,
                    " ",
                    pdbid,
                    " ",
                    chain,
                    " : MSA sequence residue at ",
                    i,
                    " (",
                    seq[i],
                    ") != PDB residue at ",
                    pdb_resnum,
                    " (",
                    pdb_res,
                    ")",
                )
                throw(ErrorException(msg))
            end
        end
    end
    m
end

function msacolumn2pdbresidue(
    msa::AnnotatedMultipleSequenceAlignment,
    seqid::String,
    pdbid::String,
    chain::String,
    pfamid::String;
    kargs...,
)
    msacolumn2pdbresidue(msa, seqid, pdbid, chain, pfamid, downloadsifts(pdbid); kargs...)
end

function msacolumn2pdbresidue(
    msa::AnnotatedMultipleSequenceAlignment,
    seqid::String,
    pdbid::String,
    chain::String;
    kargs...,
)
    msacolumn2pdbresidue(
        msa,
        seqid,
        pdbid,
        chain,
        String(split(getannotfile(msa, "AC"), '.')[1]);
        kargs...,
    )
end

"""
Returns a `BitVector` where there is a `true` for each column with PDB residue.
"""
function hasresidues(
    msa::AnnotatedMultipleSequenceAlignment,
    column2residues::AbstractDict{Int,String},
)
    colmap = getcolumnmapping(msa)
    ncol = length(colmap)
    mask = falses(ncol)
    for i = 1:ncol
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
function msaresidues(
    msa::AnnotatedMultipleSequenceAlignment,
    residues::AbstractDict{String,PDBResidue},
    column2residues::AbstractDict{Int,String},
)
    colmap = getcolumnmapping(msa)
    msares = sizehint!(OrderedDict{Int,PDBResidue}(), length(colmap))
    for col in colmap
        resnum = get(column2residues, col, "")
        if resnum != ""
            if haskey(residues, resnum)
                msares[col] = residues[resnum]
            else
                @warn(
                    "MSA column $col : The residue number $resnum isn't in the residues Dict."
                )
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
function msacontacts(
    msa::AnnotatedMultipleSequenceAlignment,
    residues::AbstractDict{String,PDBResidue},
    column2residues::AbstractDict{Int,String},
    distance_limit::Float64 = 6.05,
)
    colmap = getcolumnmapping(msa)
    contacts = columnpairsmatrix(msa)
    plm = getarray(contacts)
    @inbounds @iterateupper plm false begin

        resi = get(column2residues, colmap[i], "")
        resj = get(column2residues, colmap[j], "")
        if resi != "" && resj != "" && haskey(residues, resi) && haskey(residues, resj)
            list[k] = Float64(contact(residues[resi], residues[resj], distance_limit))
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
function getcontactmasks(contact_list::Vector{T}) where {T<:AbstractFloat}
    N = length(contact_list)
    true_contacts = falses(N)
    false_contacts = falses(N)
    @inbounds for i = 1:N
        value = contact_list[i]
        if value == 1.0
            true_contacts[i] = true
        elseif value == 0.0
            false_contacts[i] = true
        end
        # If value is NaN, It keeps the false value
    end
    true_contacts, false_contacts
end

function getcontactmasks(plm::PairwiseListMatrix{T,false,VT}) where {T<:AbstractFloat,VT}
    getcontactmasks(getlist(plm))
end

function getcontactmasks(
    nplm::NamedArray{T,2,PairwiseListMatrix{T,false,TV},DN},
) where {T,TV,DN}
    getcontactmasks(getarray(nplm))
end
