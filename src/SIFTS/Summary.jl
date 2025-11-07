# This file has the utilities needed to download and parse the CSV summary files from SIFTS
# That contains chain-level mappings between PDB and other databases.

struct SIFTSCSV <: FileFormat end

# Download Summary Files
# =======================

const _CSV_URL = "https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/"

@inline _summary_url(::Type{dbUniProt}) = _CSV_URL * "pdb_chain_uniprot.csv.gz"

@inline _summary_url(::Type{dbPfam}) = _CSV_URL * "pdb_chain_pfam.csv.gz"

@inline _summary_url(::Type{dbInterPro}) = _CSV_URL * "pdb_chain_interpro.csv.gz"

@inline _summary_url(::Type{dbSCOP}) = _CSV_URL * "pdb_chain_scop_uniprot.csv.gz"

@inline _summary_url(::Type{dbSCOP2}) = _CSV_URL * "pdb_chain_scop2_uniprot.csv.gz"

@inline _summary_url(::Type{dbSCOP2B}) = _CSV_URL * "pdb_chain_scop2b_sf_uniprot.csv.gz"

@inline _summary_url(::Type{dbCATH}) = _CSV_URL * "pdb_chain_cath_uniprot.csv.gz"

@inline _summary_url(::Type{dbEnsembl}) = _CSV_URL * "pdb_chain_ensembl.csv.gz"

@inline function _summary_url(::Type{T}) where T<:DataBase
    throw(ErrorException("No SIFTS summary file url defined for database type $T"))
end

function downloadsifts(
    database::Type{T};
    filename::Union{AbstractString,Nothing} = nothing,
) where {T<:DataBase}
    url = _summary_url(database)
    path = isnothing(filename) ? basename(url) : filename
    download_file(url, path)
end

# Parsing Summary Files
# =====================

function Utils.parse_file(io::IO, ::Type{SIFTSCSV})
    data, header = readdlm(io, ',', String, comments = true, header = true)
    colnames = Symbol.(vec(header))
    return (colnames = colnames, table = data)
end

