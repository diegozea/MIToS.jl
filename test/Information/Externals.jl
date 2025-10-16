@testset "_create_script" begin
    base = tempname()
    script = base * ".jl"
    msa_file = base * ".fasta"
    jdl_file = base * ".jls"
    MIToS.Information._create_script(
        script,
        msa_file,
        jdl_file;
        min_separation = 2,
        mode = :fast,
    )
    content = read(script, String)
    @test occursin("min_separation=2", content)
    @test occursin("mode=:fast", content)
    @test occursin("GaussDCA.gDCA(\"$msa_file\"", content)
    @test_throws AssertionError MIToS.Information._create_script(
        script,
        msa_file,
        jdl_file;
        bad = "x",
    )
    rm(script)
end

if VERSION >= v"1.5.0"
    local gaussdca_installed = false
    try
        using Pkg
        Pkg.add(
            PackageSpec(
                url = "https://github.com/carlobaldassi/GaussDCA.jl",
                rev = "master",
            ),
        )
        gaussdca_installed = true
    catch err
        @warn "GaussDCA.jl not gaussdca_installed: $err"
    end

    if gaussdca_installed
        msa = map(Residue, rand(1:21, 100, 20))
        dca = gaussdca(msa, min_separation = 2)

        @test isnan(dca[1, 1])
        @test isnan(dca[1, 2])
        @test !isnan(dca[1, 3])
    end
end
