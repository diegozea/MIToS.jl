# This file has the utilities needed to download and parse the CSV summary files from SIFTS
# That contains chain-level mappings between PDB and other databases.

# NOTE: I am using dbSCOP2 for the example as it is currently (07/11/2025) the smallest
# file (around 849KB) at https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/
"""
    SIFTSCSV

A `FileFormat` subtype for the chain-level CSV summary tables described in the
[PDBe SIFTS Quick Access guide](https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html).
Use it with the `read_file` function to load gzipped CSV files downloaded from SIFTS via
[`downloadsifts`](@ref). For example, to download and read the summary file for SCOP2:

```julia
using MIToS.SIFTS
summary_path = downloadsifts(dbSCOP2)
summary = read_file(summary_path, SIFTSCSV)
```
"""
struct SIFTSCSV <: FileFormat end

# Download Summary Files
# =======================

const _CSV_URL = "https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/"

@inline _summary_name(::Type{dbUniProt}) = "pdb_chain_uniprot.csv.gz"
@inline _summary_name(::Type{dbPfam}) = "pdb_chain_pfam.csv.gz"
@inline _summary_name(::Type{dbInterPro}) = "pdb_chain_interpro.csv.gz"
@inline _summary_name(::Type{dbSCOP}) = "pdb_chain_scop_uniprot.csv.gz"
@inline _summary_name(::Type{dbSCOP2}) = "pdb_chain_scop2_uniprot.csv.gz"
@inline _summary_name(::Type{dbSCOP2B}) = "pdb_chain_scop2b_sf_uniprot.csv.gz"
@inline _summary_name(::Type{dbCATH}) = "pdb_chain_cath_uniprot.csv.gz"
@inline _summary_name(::Type{dbEnsembl}) = "pdb_chain_ensembl.csv.gz"

@inline _summary_url(db::Type{dbUniProt}) = _CSV_URL * _summary_name(db)
@inline _summary_url(db::Type{dbPfam}) = _CSV_URL * _summary_name(db)
@inline _summary_url(db::Type{dbInterPro}) = _CSV_URL * _summary_name(db)
@inline _summary_url(db::Type{dbSCOP}) = _CSV_URL * _summary_name(db)
@inline _summary_url(db::Type{dbSCOP2}) = _CSV_URL * _summary_name(db)
@inline _summary_url(db::Type{dbSCOP2B}) = _CSV_URL * _summary_name(db)
@inline _summary_url(db::Type{dbCATH}) = _CSV_URL * _summary_name(db)
@inline _summary_url(db::Type{dbEnsembl}) = _CSV_URL * _summary_name(db)

@inline function _summary_url(::Type{T}) where T<:DataBase
    throw(ErrorException("No SIFTS summary file url defined for database type $T"))
end

"""
    downloadsifts(database::Type{<:DataBase}; filename=nothing)

Download a SIFTS chain-level summary file (CSV, gzipped) from PDBe into the
current working directory.

The `database` argument selects which summary file to fetch; see the
[PDBe SIFTS Quick Access guide](https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html)
for details on each file. This function always downloads the gzipped CSV variant.

For example, to download the `"pdb_chain_scop2_uniprot.csv.gz"` file with the chain-level
SCOP2 mappings, use `downloadsifts(dbSCOP2)`. Then, to read and parse that file,
you can use `read_file` with the [`SIFTSCSV`](@ref) format.

If `filename` is not provided, the canonical PDBe filename is used; otherwise,
the data are saved to the specified path. The function returns the path to the
downloaded file.
"""
function downloadsifts(
    database::Type{T};
    filename::AbstractString = _summary_name(database),
) where {T<:DataBase}
    url = _summary_url(database)
    download_file(url, filename)
end

# Parsing Summary Files
# =====================

"""
    Utils.parse_file(io::IO, ::Type{SIFTSCSV})

Parse a SIFTS summary CSV file from an already-open, decompressed `io` stream.
This function **expects** that `io` yields the plain CSV text (i.e., any `.gz`
decompression has already been performed). This is automatically handled if you
call `read_file` with [`SIFTSCSV`](@ref) as the format; as `read_file` opens the file, 
handles decompression when necessary, and then calls `parse_file`.

Returns a `NamedTuple` with:
- `colnames` — a `Vector{Symbol}` of column names
- `table` — the raw `Matrix{String}` produced by `DelimitedFiles.readdlm`

This low-level representation is intended for downstream reshaping or
conversion. For example, if you have the output of this function stored in
the variable `summary`, you can convert it to a `DataFrame` as follows:

```julia
using DataFrames
summary_df = DataFrame(summary.table, summary.colnames)
```
"""
function Utils.parse_file(io::IO, ::Type{SIFTSCSV})
    data, header = readdlm(io, ',', String, comments = true, header = true, quotes=false)
    # quotes=false is needed to parse the taxonomy file correctly
    colnames = Symbol.(vec(header))
    return (colnames = colnames, table = data)
end

