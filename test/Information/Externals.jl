if VERSION >= v"1.5.0"
    installed = false
    try
        using Pkg
        Pkg.add(PackageSpec(url="https://github.com/carlobaldassi/GaussDCA.jl", rev="master"))
        installed = true
    catch err
        @warn "GaussDCA.jl not installed: $err"
    end

    if installed
        msa = map(Residue, rand(1:21,100,20))
        dca = gaussdca(msa, min_separation=2)

        @test  isnan(dca[1,1])
        @test  isnan(dca[1,2])
        @test !isnan(dca[1,3])
    end
end