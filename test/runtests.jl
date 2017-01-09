using Base.Test

using MIToS: Utils, MSA, Information, PDB

using NamedArrays           # array
using StatsBase             # WeightVec
using PairwiseListMatrices  # getlist

# using MIToS.Information
# using MIToS.PDB
# using MIToS.SIFTS
# using MIToS.Pfam

const DATA = joinpath(pwd(), "data")

# # Utils
# @testset "Utils" begin
#     include(joinpath("Utils", "GeneralUtils.jl"))
# end
#
# # MSA
# @testset "MSA" begin
#     include(joinpath("MSA", "Residues.jl"))
#     include(joinpath("MSA", "Alphabet.jl"))
#     include(joinpath("MSA", "ThreeLetters.jl"))
#     include(joinpath("MSA", "Annotations.jl"))
#     include(joinpath("MSA", "MultipleSequenceAlignment.jl"))
#     include(joinpath("MSA", "GeneralParserMethods.jl"))
#     include(joinpath("MSA", "IO.jl"))
#     include(joinpath("MSA", "General.jl"))
#     include(joinpath("MSA", "MSAEditing.jl"))
#     include(joinpath("MSA", "MSAStats.jl"))
#     include(joinpath("MSA", "Shuffle.jl"))
#     include(joinpath("MSA", "Identity.jl"))
#     include(joinpath("MSA", "Hobohm.jl"))
#     include(joinpath("MSA", "MSAAnnotations.jl"))
# end
#
# # Information
# @testset "Information" begin
#     include(joinpath("Information", "ContingencyTables.jl"))
#     include(joinpath("Information", "Counters.jl"))
#     include(joinpath("Information", "InformationMeasures.jl"))
#     include(joinpath("Information", "Iterations.jl"))
#     include(joinpath("Information", "CorrectedMutualInformation.jl"))
#     include(joinpath("Information", "Gaps.jl"))
# end

# PDB
@testset "PDB" begin
    include(joinpath("PDB", "PDB.jl"))
    include(joinpath("PDB", "Contacts.jl"))
    include(joinpath("PDB", "Kabsch.jl"))
end

# # MSA
# include("indexedarrays.jl") # TO DO: Agregar si es utilizado luego
# # PDB
# include("kabsch.jl")
# #  SIFTS
# include("sifts.jl")
# # Pfam
# include("pfam.jl")
# # Scripts
# include("scripts.jl")

print("""

----- =D -----

""")
