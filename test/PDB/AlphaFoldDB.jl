@testset "AlphaFoldDB Tests" begin
    @testset "query_alphafolddb tests" begin
        @testset "Valid UniProt ID" begin
            structure_info = query_alphafolddb("A0A0C5B5G6")
            @test structure_info[1]["uniprotAccession"] == "A0A0C5B5G6"
            @test structure_info[1]["uniprotId"] == "MOTSC_HUMAN"
            @test match(r"^https.+\.pdb$", structure_info[1]["pdbUrl"]) !== nothing
        end

        # @testset "Invalid UniProt ID" begin
        #     @test_throws HTTP.Exceptions.StatusError query_alphafolddb("INVALID_ID")
        # end
    end

    @testset "download_alphafold_structure tests" begin
        @testset "Download PDB format" begin
            outfile = download_alphafold_structure("A0A0C5B5G6", format=PDBFile)
            if isfile(outfile)
                try
                    res = read_file(outfile, PDBFile)
                    @test length(res) == 16
                finally
                    rm(outfile)
                end
            end
        end

        @testset "Unsupported format" begin
            struct UnsupportedFormat <: FileFormat end
            @test_throws ArgumentError download_alphafold_structure("A0A0C5B5G6", format=UnsupportedFormat)
        end
    end

    @testset "_extract_filename_from_url tests" begin
        @testset "Extract filename" begin
            filename = MIToS.PDB._extract_filename_from_url("https://example.com/structure.pdb")
            @test filename == "structure.pdb"
        end
    end
end
