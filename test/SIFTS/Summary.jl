@testset "SIFTS Summary" begin

    implemented_databases = Set([
        dbUniProt,
        dbPfam,
        dbInterPro,
        dbSCOP,
        dbSCOP2,
        dbSCOP2B,
        dbCATH,
        dbEnsembl,
    ])
    
    unimplemented_databases = setdiff(Set(subtypes(DataBase)), implemented_databases)
    # i.e., all the DataBase subtypes that are not included in the implemented_databases

    @testset "URL generation" begin
        for db in unimplemented_databases
            @test_throws ErrorException MIToS.SIFTS._summary_url(db)
        end
        for db in implemented_databases
            url = MIToS.SIFTS._summary_url(db)
            @test startswith(url, "https://")
            @test endswith(url, ".csv.gz")
        end
    end

    # We will test using only dbSCOP2 as it is currently (07/11/2025) the smallest
    # file (around 849KB) at https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/
    @testset "Downloading and parsing" begin
        mktempdir() do dir
            cd(dir) do
                # Test Download
                downloadsifts(dbSCOP2)
                @test isfile("pdb_chain_scop2_uniprot.csv.gz")
                # Test Parsing
                data = read_file("pdb_chain_scop2_uniprot.csv.gz", SIFTSCSV)
                @test isa(data, NamedTuple)
                @test haskey(data, :colnames)
                @test haskey(data, :table)
                @test isa(data.colnames, Vector{Symbol})
                @test isa(data.table, Matrix{String})
                @test length(data.colnames) >= 3 # Ensure at least 3 columns are present
                @test size(data.table, 2) == length(data.colnames) # Ensure columns match
                @test size(data.table, 1) >= 10000 # Ensure at least 10,000 rows are present
                # NOTE: We do not test exact numbers as the files may be updated over time
            end
        end
    end
end