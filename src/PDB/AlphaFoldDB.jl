using Downloads
using JSON3
using HTTP
using MIToS.PDB
using MIToS.Utils

"""
    query_alphafolddb(uniprot_id::String)

This function queries the AlphaFoldDB API to retrieve structure information for 
a given `uniprot_id`. This function returns the structure information as a
JSON object.
"""
function query_alphafolddb(uniprot_id::String)
    # Construct the URL for the AlphaFoldDB API request
    url = "https://alphafold.ebi.ac.uk/api/prediction/$uniprot_id"

    response = HTTP.request("GET", url)

    if response.status == 200
        # Read the JSON response body
        JSON3.read(String(response.body))
    else
        throw(ErrorException(
            "Error fetching UniProt ID $uniprot_id from AlphaFoldDB. Status: $(response.status)"))
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

# Example of using the function with a given UniProt ID
uniprot_id = "Q5VSL9"
download_alphafold_structure(uniprot_id)


# ---

#=
f(x) = x
f(x::Int) = "$x is an integer"
f(x::AbstractFloat) = "$x is a float"
f(x::Int, y::Int) = "$x and $y are integers"
f(x::Int, y::AbstractFloat) = "$x is an integer and $y is a float"
f(x::Int, y::String) = "$x is an integer and $y is a string"
f(x, y) = nothing

g(x::T) where T <: AbstractFloat = "$x is a $T"
=#
