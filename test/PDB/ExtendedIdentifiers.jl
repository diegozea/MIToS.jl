@testset "Extended PDB identifiers" begin
    @testset "Download paths and errors" begin
        for (code, name) in [
            ("1AbC", "1ABC"),
            ("PDB_00001ABC", "pdb_00001abc"),
            ("pDb_10021AbC", "pdb_10021abc"),
        ]
            for (format, extension) in [(MMCIFFile, ".cif.gz"), (PDBML, ".xml.gz")]
                @test PDB._pdb_download_name(code, format) == name * extension
                for baseurl in ["https://example.org/files", "https://example.org/files/"]
                    @test PDB._pdb_download_url(code, format, baseurl) ==
                          "https://example.org/files/" * name * extension
                end
            end
        end
        @test PDB._pdb_baseurl("1abc") == "https://files.rcsb.org/download/"
        @test PDB._pdb_baseurl("pdb_00001abc") == "https://files-beta.wwpdb.org/download/"
        @test PDB._pdb_download_name("PDB_00001ABC", PDBFile) == "pdb_00001abc.pdb.gz"
        @test PDB._graphql_query("PDB_00001ABC") == PDB._graphql_query("1abc")

        mktempdir() do dir
            path = joinpath(dir, "existing.gz")
            write(path, "keep this file")
            for code in ["pdb_10021abc", "pdb_00011a7y", "pdb_0000abcd"]
                @test_throws ArgumentError downloadpdb(
                    code,
                    format = PDBFile,
                    filename = path,
                )
                @test_throws ArgumentError getpdbdescription(code)
                @test_throws ArgumentError downloadpdbheader(code, filename = path)
                @test read(path, String) == "keep this file"
            end
            for code in ["1ab_", "pdb_00001abc\n", "../1abc"]
                @test_throws ErrorException downloadpdb(code, filename = path)
                @test_throws ErrorException downloadpdbheader(code, filename = path)
                @test read(path, String) == "keep this file"
            end
        end
    end

    @testset "Official wwPDB extended-only example" begin
        path = joinpath(DATA, "pdb_00011a7y.cif.gz")
        residues = read_file(path, MMCIFFile)
        @test length(residues) == 41
        @test sum(length(res.atoms) for res in residues) == 314
        @test count(res -> res.id.name == "X2AVD", residues) == 6
        auth_residues = read_file(path, MMCIFFile, label = false)
        @test length(auth_residues) == 41
        @test count(res -> res.id.name == "X2AVD", auth_residues) == 6
        # The fixture contains a real extended data block and entry ID, not just a new filename.
        open(path) do io
            stream = Utils.GzipDecompressorStream(io)
            dict = BioStructures.MMCIFDict(stream)
            @test dict["data_"] == ["pdb_00011a7y"]
            @test dict["_entry.id"] == ["pdb_00011a7y"]
            close(stream)
        end
    end

    @testset "Beta archive downloads" begin
        mktempdir() do dir
            cd(dir) do
                for (format, extension) in
                    [(MMCIFFile, ".cif.gz"), (PDBML, ".xml.gz"), (PDBFile, ".pdb.gz")]
                    filename = downloadpdb("PDB_00002VQC", format = format)
                    @test filename == "pdb_00002vqc" * extension
                    residues = read_file(filename, format)
                    @test findfirst(res -> res.id.number == "4", residues) == 1
                    @test findfirst(res -> res.id.number == "73", residues) == 70
                end
                # A custom mirror and destination still work; RCSB accepts extended mmCIF IDs.
                filename = downloadpdb(
                    "pdb_00002vqc",
                    filename = "custom.cif",
                    baseurl = "https://files.rcsb.org/download",
                )
                @test filename == "custom.cif.gz"
                @test !isempty(read_file(filename, MMCIFFile))
            end
        end
        @test getpdbdescription("pDb_00004hHb")["rcsb_entry_info"]["resolution_combined"][1] ==
              1.74
    end
end
