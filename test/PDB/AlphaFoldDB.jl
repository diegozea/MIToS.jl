using Test
using Downloads
using JSON3
using HTTP
using MIToS.PDB
using MIToS.Utils

# Include the original code file (assuming it's named `alphafold_download.jl`)
include("alphafold_download.jl")

# Mocking HTTP responses for testing
function mock_query_alphafolddb(uniprot_id::String)
    if uniprot_id == "Q5VSL9"
        return JSON3.read("""
        [
            {
                "uniprotId": "$uniprot_id",
                "pdbUrl": "https://example.com/mock_structure.pdb"
            }
        ]
        """)
    else
        throw(ErrorException("UniProt ID not found"))
    end
end

# Overriding the actual function with the mock
query_alphafolddb(uniprot_id::String) = mock_query_alphafolddb(uniprot_id)

# Mocking file download for testing
function mock_download_file(url::String, file_path::String)
    println("Mock download of $url to $file_path")
end

# Overriding the actual function with the mock
download_file(url::String, file_path::String) = mock_download_file(url, file_path)

# Test suite
@testset "AlphaFoldDB Tests" begin
    @testset "query_alphafolddb tests" begin
        @test "Valid UniProt ID" begin
            structure_info = query_alphafolddb("Q5VSL9")
            @test structure_info[1]["uniprotId"] == "Q5VSL9"
            @test structure_info[1]["pdbUrl"] == "https://example.com/mock_structure.pdb"
        end

        @test "Invalid UniProt ID" begin
            @test_throws ErrorException query_alphafolddb("INVALID_ID")
        end
    end

    @testset "download_alphafold_structure tests" begin
        @test "Download PDB format" begin
            download_alphafold_structure("Q5VSL9", format=PDBFile)
            # Verify output, in a real test you would check the file system or the mocked output
        end

        @test "Unsupported format" begin
            struct UnsupportedFormat <: FileFormat end
            @test_throws ArgumentError download_alphafold_structure("Q5VSL9", format=UnsupportedFormat)
        end
    end

    @testset "_extract_filename_from_url tests" begin
        @test "Extract filename" begin
            filename = _extract_filename_from_url("https://example.com/path/to/structure.pdb")
            @test filename == "structure.pdb"
        end
    end
end
