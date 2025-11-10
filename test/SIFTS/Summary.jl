@testset "SIFTS Summary" begin

    implemented_databases =
        Set([dbUniProt, dbPfam, dbInterPro, dbSCOP, dbSCOP2, dbSCOP2B, dbCATH, dbEnsembl])

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

    @testset "Downloading and parsing" begin
        mktempdir() do dir
            cd(dir) do
                # This is time consuming, but we need to ensure that downloadsifts and
                # and read_file with SIFTSCSV work as expected for all supported databases.
                for db in implemented_databases
                    # Test Download
                    csv_file = downloadsifts(db)
                    # Test Parsing
                    data = read_file(csv_file, SIFTSCSV)
                    @test isa(data, NamedTuple)
                    @test haskey(data, :colnames)
                    @test haskey(data, :table)
                    @test isa(data.colnames, Vector{Symbol})
                    @test isa(data.table, Matrix{String})
                    @test length(data.colnames) >= 3 # Ensure at least 3 columns are present
                    @test size(data.table, 2) == length(data.colnames) # Ensure columns match
                    @test size(data.table, 1) >= 5 # Ensure at least 5 rows are present
                end
            end
        end
    end
end