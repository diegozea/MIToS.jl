@testset "Scripts" begin

    # julia bin
    julia = joinpath(Base.Sys.BINDIR, "julia")
    # ../../
    mitos_folder = splitdir(splitdir(dirname(@__FILE__))[1])[1]
    # current project
    project = Base.active_project()

    @testset "Distances.jl" begin
        path_script = joinpath(mitos_folder, "scripts", "Distances.jl")
        path_file = joinpath(mitos_folder, "test", "data", "small.pdb")
        out_intra =
            read(`$julia --project=$project $path_script $path_file -o STDOUT`, String)
        out_inter = read(
            `$julia --project=$project $path_script $path_file --inter -o STDOUT`,
            String,
        )

        intra = Float64[
            parse(Float64, split(line, ',')[end]) for line in split(out_intra, '\n') if
            !startswith(line, '#') && occursin(r"\d\.\d+$", line)
        ]
        @test sum(intra) == 2.0 * 3.7836265671971385
        @test length(collect((m.match for m in eachmatch(r"\n1,A,", out_intra)))) == 1
        @test length(collect((m.match for m in eachmatch(r"\n1,B,", out_intra)))) == 1

        inter = Float64[
            parse(Float64, split(line, ',')[end]) for line in split(out_inter, '\n') if
            !startswith(line, '#') && occursin(r"\d\.\d+$", line)
        ]
        @test sum(inter) == 4.0 * 3.7836265671971385
        @test length(collect((m.match for m in eachmatch(r"\n1,A,", out_inter)))) == 5
        @test length(collect((m.match for m in eachmatch(r"\n1,B,", out_inter)))) == 1
    end

    @testset "AlignedColumns.jl" begin
        path_script = joinpath(mitos_folder, "scripts", "AlignedColumns.jl")
        path_file = joinpath(mitos_folder, "test", "data", "PF09645_full.stockholm")
        output = read(`$julia --project=$project $path_script $path_file -o STDOUT`, String)

        @test read_file(
            path_file,
            Stockholm,
            generatemapping = true,
            useidcoordinates = true,
            deletefullgaps = true,
        ) == parse_file(output, Stockholm)
    end

    @testset "BLMI.jl" begin
        path_script = joinpath(mitos_folder, "scripts", "BLMI.jl")
        path_file = joinpath(mitos_folder, "test", "data", "simple.fasta")
        output = read(
            `$julia --project=$project $path_script $path_file -o STDOUT --samples 0 --format FASTA --apc`,
            String,
        )

        @test occursin(r"1,2,0.0,0.17", output)
    end

    @testset "Buslje09.jl" begin
        path_script = joinpath(mitos_folder, "scripts", "Buslje09.jl")
        path_file = joinpath(mitos_folder, "test", "data", "simple.fasta")
        output = read(
            `$julia --project=$project $path_script $path_file -o STDOUT --samples 0 --format FASTA --apc`,
            String,
        )

        @test occursin(r"1,2,0.0,0.13", output)
    end

    @testset "Conservation.jl" begin
        path_script = joinpath(mitos_folder, "scripts", "Conservation.jl")
        path_file = joinpath(mitos_folder, "test", "data", "simple.fasta")
        output = read(
            `$julia --project=$project $path_script $path_file -o STDOUT -c --format FASTA`,
            String,
        )

        @test occursin(r"1,0.6931471805599453,2.0901394960274127", output)
        @test occursin(r"2,0.6931471805599453,2.0901394960274127", output)
    end

    @testset "DownloadPDB.jl" begin
        # Helper functions
        _get_files() = filter!(f -> startswith(f, "3NIR"), readdir())
        function _delete_files()
            for file in _get_files()
                rm(file)
            end
        end
        # Test
        _delete_files()
        path_script = joinpath(mitos_folder, "scripts", "DownloadPDB.jl")
        output = read(`$julia --project=$project $path_script --code 3NIR`, String)
        files = _get_files()
        @test length(files) == 1
        _delete_files()
    end

    @testset "MSADescription.jl" begin
        path_script = joinpath(mitos_folder, "scripts", "MSADescription.jl")
        path_file = joinpath(mitos_folder, "test", "data", "PF09645_full.stockholm")
        output = read(`$julia --project=$project $path_script $path_file -o STDOUT`, String)

        @test occursin(r"clusters,number,,4", output)
        @test occursin(r"gapfraction,quantile,0.75,0.25", output)
        @test occursin(r"coverage,mean,,0.80", output)
    end

    @testset "PairwiseGapPercentage.jl" begin
        path_script = joinpath(mitos_folder, "scripts", "PairwiseGapPercentage.jl")
        path_file = joinpath(mitos_folder, "test", "data", "PF09645_full.stockholm")
        output = read(`$julia --project=$project $path_script $path_file -o STDOUT`, String)

        @test occursin(r"115,115,75.0,75.0", output)
        @test occursin(r"40,70,0.0,0.0", output)
        @test occursin(r"115,115,75.0,75.0", output)
    end

    @testset "PercentIdentity.jl" begin
        path_script = joinpath(mitos_folder, "scripts", "PercentIdentity.jl")
        path_file = joinpath(mitos_folder, "test", "data", "PF09645_full.stockholm")
        output = read(`$julia --project=$project $path_script $path_file -o STDOUT`, String)

        @test occursin(
            r"110,4,29.5,15.21,14.13,19.75[0-9]*,26.3[0-9]+,33.78[0-9]*,56.38",
            output,
        )
    end
end
