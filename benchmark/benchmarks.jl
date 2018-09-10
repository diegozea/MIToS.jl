using PkgBenchmark
using Random
using MIToS.Utils
using MIToS.MSA

@benchgroup "Utils" begin
    include(joinpath("Utils", "GeneralUtils.jl"))
end

@benchgroup "MSA" begin
    include(joinpath("MSA", "Residues.jl"))
    include(joinpath("MSA", "Annotations.jl"))
end
