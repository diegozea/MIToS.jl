using BenchmarkTools
using Random
using MIToS.Utils
using MIToS.MSA
using MIToS.Information

const SUITE = BenchmarkGroup()

include("Utils/GeneralUtils.jl")
include("MSA/Residues.jl")
include("MSA/Annotations.jl")
include("Information/CorrectedMutualInformation.jl")