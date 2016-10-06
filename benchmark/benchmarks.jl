using PkgBenchmark
using MIToS: Utils

@benchgroup "Utils" begin
    include(joinpath("Utils", "GeneralUtils.jl"))
end
