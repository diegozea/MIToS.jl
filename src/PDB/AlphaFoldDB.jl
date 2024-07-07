"""
    query_alphafolddb(uniprot_accession::String)

This function queries the AlphaFoldDB API to retrieve structure information for
a given `uniprot_accession`, e.g. `"P00520"`. This function returns the structure
information as a `JSON3.Object`.
"""
function query_alphafolddb(uniprot_accession::String)
    # Construct the URL for the AlphaFoldDB API request
    url = "https://alphafold.ebi.ac.uk/api/prediction/$uniprot_accession"

    body = IOBuffer()
    response = Downloads.request(url, method = "GET", output = body)

    if response.status == 200
        # Use only to get the unique EntrySummary object in the Root list
        only(JSON3.read(String(take!(body))))
    else
        error_type = response.status == 422 ? "Validation Error" : "Error"
        throw(
            ErrorException(
                "$error_type ($(response.status)) fetching UniProt Accession $uniprot_accession from AlphaFoldDB.",
            ),
        )
    end
end

# This function extracts the filename from a given URL.
function _extract_filename_from_url(url::String)
    return split(url, "/")[end]
end

# Function to download the PDB or CIF file based on the UniProt Accession
"""
    download_alphafold_structure(uniprot_accession::String; format::Type{T}=MMCIFFile) where T<:FileFormat

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

    # Initialize the model URL based on the requested format
    if format === PDBFile
        model_url = structure_info["pdbUrl"]
    elseif format === MMCIFFile
        model_url = structure_info["cifUrl"]
    else
        throw(ArgumentError("Unsupported format: $format"))
    end

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
