"""
    query_alphafolddb(uniprot_accession::String)

This function queries the AlphaFoldDB API to retrieve structure information for
a given `uniprot_accession`, e.g. `"P00520"`. This function returns the structure
information as a `JSON3.Object`. The download is performed via `download_file`,
which already retries the request using an exponential backoff strategy.
"""
function query_alphafolddb(uniprot_accession::AbstractString)
    # Construct the URL for the AlphaFoldDB API request
    url = "https://alphafold.ebi.ac.uk/api/prediction/$uniprot_accession"

    filepath = try
        download_file(url)
    catch err
        @error "Error fetching UniProt Accession $uniprot_accession from AlphaFoldDB."
        rethrow(err)
    end
    try
        body = read(filepath, String)
        return only(JSON3.read(body))
    catch err
        @error "Unexpected AlphaFoldDB response for $uniprot_accession."
        rethrow(err)
    finally
        isfile(filepath) && rm(filepath)
    end
end

# This function extracts the filename from a given URL.
function _extract_filename_from_url(url::AbstractString)
    return split(url, "/")[end]
end

# Function to download the PDB or CIF file based on the UniProt Accession
"""
    download_alphafold_structure(uniprot_accession::AbstractString; format::Type{T}=MMCIFFile) where T<:FileFormat

This function downloads the structure file (PDB or mmCIF) for a given UniProt Accession
from AlphaFoldDB. The `uniprot_accession` parameter specifies the UniProt Accession of the
protein, e.g. `"P00520"`. The `format` parameter specifies the file format to download,
with the default being mmCIF, i.e. `MMCIFFile`. You can set `format` to `PDBFile` if you
want to download a PDB file.
"""
function download_alphafold_structure(
    uniprot_accession::String;
    format::Type{T} = MMCIFFile,
) where {T<:FileFormat}

    structure_info = query_alphafolddb(uniprot_accession)
    @assert format === PDBFile || format === MMCIFFile "Unsupported format: $format"
    # Initialize the model URL based on the requested format
    model_url = format === PDBFile ? structure_info["pdbUrl"] : structure_info["cifUrl"]

    file_name = _extract_filename_from_url(model_url)

    try
        download_file(model_url, file_name)
    catch
        throw(
            ErrorException(
                "Error downloading AlphaFold model for UniProt Accession $uniprot_accession",
            ),
        )
    end
end
