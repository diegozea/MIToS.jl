@testset "Scripts" begin

    # julia bin
    julia = joinpath(Base.JULIA_HOME, "julia")
    # ../../
    mitos_folder = splitdir(splitdir(dirname(@__FILE__))[1])[1]

    @testset "Distances.jl" begin
        path_script = joinpath(mitos_folder, "scripts", "Distances.jl")
        path_file   = joinpath(mitos_folder, "test", "data", "small.pdb")
        out_intra = readstring(`$julia $path_script $path_file -o STDOUT`)
        out_inter = readstring(`$julia $path_script $path_file --inter -o STDOUT`)

        intra = Float64[parse(Float64,split(line,',')[end]) for line in
            split(out_intra,'\n') if !startswith(line,'#') && ismatch(r"\d\.\d+$",line)]
        @test sum(intra) == 2.0*3.7836265671971385
        @test length(matchall(r"\n1,A,", out_intra)) == 1
        @test length(matchall(r"\n1,B,", out_intra)) == 1

        inter = Float64[parse(Float64,split(line,',')[end]) for line in
            split(out_inter,'\n') if !startswith(line,'#') && ismatch(r"\d\.\d+$",line)]
        @test sum(inter) == 4.0*3.7836265671971385
        @test length(matchall(r"\n1,A,", out_inter)) == 5
        @test length(matchall(r"\n1,B,", out_inter)) == 1
    end

    @testset "AlignedColumns.jl" begin
        path_script = joinpath(mitos_folder, "scripts", "AlignedColumns.jl")
        path_file   = joinpath(mitos_folder, "test", "data", "PF09645_full.stockholm")
        output = readstring(`$julia $path_script $path_file -o STDOUT`)

        @test read(path_file, Stockholm, generatemapping=true, useidcoordinates=true,
            deletefullgaps=true) == parse(output, Stockholm)
    end

    @testset "BLMI.jl" begin
        path_script = joinpath(mitos_folder, "scripts", "BLMI.jl")
        path_file   = joinpath(mitos_folder, "test", "data", "simple.fasta")
        output = readstring(`$julia $path_script $path_file -o STDOUT --samples 0 --format FASTA --apc`)

        @test ismatch(r"1,2,0,0.17",output)
    end

    @testset "Buslje09.jl" begin
        path_script = joinpath(mitos_folder, "scripts", "Buslje09.jl")
        path_file   = joinpath(mitos_folder, "test", "data", "simple.fasta")
        output = readstring(`$julia $path_script $path_file -o STDOUT --samples 0 --format FASTA --apc`)

        @test ismatch(r"1,2,0,0.13",output)
    end

    @testset "Conservation.jl" begin
        path_script = joinpath(mitos_folder, "scripts", "Conservation.jl")
        path_file   = joinpath(mitos_folder, "test", "data", "simple.fasta")
        output = readstring(`$julia $path_script $path_file -o STDOUT -c --format FASTA`)

        @test ismatch(r"1,0.6931471805599453,2.0901394960274127",output)
        @test ismatch(r"2,0.6931471805599453,2.0901394960274127",output)
    end

    @testset "DownloadPDB.jl" begin
        # Helper functions
        _get_files() = filter!(f -> startswith(f, "3NIR"),readdir())
        function _delete_files()
            for file in _get_files()
                rm(file)
            end
        end
        # Test
        _delete_files()
        path_script = joinpath(mitos_folder, "scripts", "DownloadPDB.jl")
        output = readstring(`$julia $path_script --code 3NIR`)
        files = _get_files()
        @test length(files) == 1
        _delete_files()
    end

    @testset "MSADescription.jl" begin
        path_script = joinpath(mitos_folder, "scripts", "MSADescription.jl")
        path_file   = joinpath(mitos_folder, "test", "data", "PF09645_full.stockholm")
        output = readstring(`$julia $path_script $path_file -o STDOUT`)

        @test ismatch(r"clusters,number,,4", output)
        @test ismatch(r"gapfraction,quantile,0.75,0.25", output)
        @test ismatch(r"coverage,mean,,0.80", output)
    end

    @testset "PairwiseGapPercentage.jl" begin
        path_script = joinpath(mitos_folder, "scripts", "PairwiseGapPercentage.jl")
        path_file   = joinpath(mitos_folder, "test", "data", "PF09645_full.stockholm")
        output = readstring(`$julia $path_script $path_file -o STDOUT`)

        @test ismatch(r"115,115,75,75", output)
        @test ismatch(r"40,70,0,0", output)
        @test ismatch(r"115,115,75,75", output)
    end

    @testset "PercentIdentity.jl" begin
        path_script = joinpath(mitos_folder, "scripts", "PercentIdentity.jl")
        path_file   = joinpath(mitos_folder, "test", "data", "PF09645_full.stockholm")
        output = readstring(`$julia $path_script $path_file -o STDOUT`)

        @test ismatch(r"110,4,29.498697,15.208886,14.13,19.75,26.31,33.78,56.38", output)
    end
end
