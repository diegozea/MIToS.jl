using BenchmarkTools
using Random
using MIToS.Utils
using MIToS.MSA
using MIToS.Information
using MIToS.PDB

const SUITE = BenchmarkGroup()

include("Utils/GeneralUtils.jl")
include("MSA/Residues.jl")
include("MSA/Annotations.jl")
include("MSA/Read.jl")
include("MSA/Write.jl")
include("MSA/Identity.jl")
include("MSA/Clustering.jl")
include("Information/CorrectedMutualInformation.jl")
include("Information/Counters.jl")
include("Information/Entropy.jl")
include("Information/HighLevel.jl")
include("Information/MIp.jl")
include("PDB/Count.jl")
include("PDB/Distance.jl")
