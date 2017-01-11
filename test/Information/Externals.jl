"Function from MLPlots.jl, written by Tom Breloff."
function _is_installed(name::String)
    try
        Pkg.installed(name) === nothing ? false : true
    catch
        false
    end
end

if _is_installed("GaussDCA")

    msa = map(Residue, rand(1:21,100,20))
    dca = gaussdca(msa, min_separation=2)

    @test  isnan(dca[1,1])
    @test  isnan(dca[1,2])
    @test !isnan(dca[1,3])
end
