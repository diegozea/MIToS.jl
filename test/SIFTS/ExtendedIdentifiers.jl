@testset "Extended SIFTS identifiers" begin
    for source in ["ftp", "https"]
        @test SIFTS._sifts_url("PDB_00002VQC", source) == SIFTS._sifts_url("2vqc", source)
        @test_throws ArgumentError SIFTS._sifts_url("pdb_10021abc", source)
    end
    @test SIFTS._sifts_url("PDB_00002VQC", "ftp") ==
          "https://ftp.ebi.ac.uk/pub/databases/msd/sifts/split_xml/vq/2vqc.xml.gz"
    @test_throws ArgumentError downloadsifts("pdb_00011a7y")

    legacy_file = joinpath(DATA, "2vqc.xml.gz")
    reference = siftsmapping(legacy_file, dbPDB, "2vqc", dbUniProt, "P20220", chain = "A")
    @test !isempty(reference)
    for code in ["2VQC", "pdb_00002vqc", "PDB_00002VQC"]
        @test siftsmapping(legacy_file, dbPDB, code, dbUniProt, "P20220", chain = "A") ==
              reference
        @test siftsmapping(legacy_file, dbUniProt, "P20220", dbPDB, code, chain = "A") ==
              siftsmapping(legacy_file, dbUniProt, "P20220", dbPDB, "2vqc", chain = "A")
    end
    # Matching of other database accessions stays case sensitive.
    @test isempty(siftsmapping(legacy_file, dbPDB, "pdb_00002vqc", dbUniProt, "p20220"))

    xml = open(legacy_file) do io
        stream = Utils.GzipDecompressorStream(io)
        try
            read(stream, String)
        finally
            close(stream)
        end
    end
    msa = read_file(
        joinpath(DATA, "PF09645_full.stockholm"),
        Stockholm,
        generatemapping = true,
        useidcoordinates = true,
    )
    seqid = "F112_SSV1/3-112"
    columns = msacolumn2pdbresidue(msa, seqid, "2VQC", "A", "PF09645", legacy_file)
    @test msacolumn2pdbresidue(msa, seqid, "PDB_00002VQC", "A", "PF09645", legacy_file) ==
          columns
    annotations = replace(
        read(joinpath(DATA, "PF09645_full.stockholm"), String),
        "PDB; 2VQC" => "PDB; pdb_10022vqc",
    )
    @test getseq2pdb(parse_file(annotations, Stockholm))[seqid] == [("pdb_10022vqc", "A")]
    summary =
        parse_file(IOBuffer("PDB,CHAIN,SP_PRIMARY\npdb_10022vqc,A,P20220\n"), SIFTSCSV)
    @test summary.table[1, 1] == "pdb_10022vqc"

    mktempdir() do dir
        path = joinpath(dir, "extended.xml")
        for code in ["PDB_00002VQC", "pdb_10022vqc"]
            # Synthetic extended SIFTS XML: exercise local parsing independently of upstream availability.
            write(path, replace(xml, "2vqc" => code))
            residues = read_file(path, SIFTSXML)
            @test all(res.PDB.id == code for res in residues if !ismissing(res.PDB))
            @test siftsmapping(path, dbPDB, code, dbUniProt, "P20220", chain = "A") ==
                  reference
            @test msacolumn2pdbresidue(msa, seqid, code, "A", "PF09645", path) == columns
            if code == "PDB_00002VQC"
                @test siftsmapping(path, dbPDB, "2vqc", dbUniProt, "P20220", chain = "A") ==
                      reference
                @test msacolumn2pdbresidue(msa, seqid, "2vqc", "A", "PF09645", path) ==
                      columns
            else
                @test isempty(siftsmapping(path, dbPDB, "2vqc", dbUniProt, "P20220"))
            end
        end
        cd(dir) do
            for source in ["ftp", "https"]
                filename = downloadsifts("PDB_00002VQC", source = source)
                @test filename == "pdb_00002vqc.xml.gz"
                mapping = siftsmapping(
                    filename,
                    dbPDB,
                    "PDB_00002VQC",
                    dbUniProt,
                    "P20220",
                    chain = "A",
                    missings = false,
                )
                @test mapping["4"] == "4"
                @test mapping == siftsmapping(
                    filename,
                    dbPDB,
                    "2vqc",
                    dbUniProt,
                    "P20220",
                    chain = "A",
                    missings = false,
                )
            end
        end
    end
end
