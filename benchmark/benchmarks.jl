using PkgBenchmark
using Random
using MIToS.Utils
using MIToS.MSA

@benchgroup "Utils" begin
    include("Utils/GeneralUtils.jl")
end

@benchgroup "MSA" begin
    include("MSA/Residues.jl")
    include("MSA/Annotations.jl")
end
