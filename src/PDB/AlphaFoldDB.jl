"""
    query_alphafolddb(uniprot_id::String)

This function queries the AlphaFoldDB API to retrieve structure information for 
a given `uniprot_id`. This function returns the structure information as a
`JSON3.Object`.
"""
function query_alphafolddb(uniprot_id::String)
    # Construct the URL for the AlphaFoldDB API request
    url = "https://alphafold.ebi.ac.uk/api/prediction/$uniprot_id"

    response = HTTP.request("GET", url)

    if response.status == 200
        JSON3.read(String(response.body))
    else
        error_type = response.status == 422 ? "Validation Error" : "Error"
        throw(ErrorException(
            "$error_type ($(response.status)) fetching UniProt ID $uniprot_id from AlphaFoldDB."))
    end
end

# This function extracts the filename from a given URL.
function _extract_filename_from_url(url::String)
    return split(url, "/")[end]
end

# Function to download the PDB or CIF file based on the UniProt ID
"""
    download_alphafold_structure(uniprot_id::String; format::Type{T}=PDBFile) where T<:FileFormat

This function downloads the structure file (PDB or CIF) for a given UniProt ID from AlphaFoldDB.
The `uniprot_id` parameter specifies the UniProt ID of the protein.
The `format` parameter specifies the file format to download, with the default being PDBFile.
"""
function download_alphafold_structure(uniprot_id::String; 
        format::Type{T}=PDBFile) where T<:FileFormat
    
    structure_info = query_alphafolddb(uniprot_id)
   
    # Initialize the model URL based on the requested format
    if format === PDBFile
        model_url = structure_info[1]["pdbUrl"]
    # elseif format === CIF
    #     model_url = structure_info[1]["cifUrl"]
    else
        throw(ArgumentError("Unsupported format: $format"))
    end
    
    file_name = _extract_filename_from_url(model_url)
    
    try
        download_file(model_url, file_name)
    catch
        throw(ErrorException("Error downloading AlphaFold model for UniProt ID $uniprot_id"))
    end
end