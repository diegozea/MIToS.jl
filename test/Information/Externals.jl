if VERSION >= v"1.5.0"
    gaussdca_installed = false
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
