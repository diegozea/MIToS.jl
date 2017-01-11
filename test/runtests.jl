using Base.Test

using MIToS: Utils, MSA, Information, PDB, SIFTS, Pfam

using NamedArrays           # array
using StatsBase             # WeightVec
using PairwiseListMatrices  # getlist
using ROCAnalysis           # AUC

const DATA = joinpath(pwd(), "data")

# Utils
@testset "Utils" begin
    include(joinpath("Utils", "GeneralUtils.jl"))
end

# MSA
@testset "MSA" begin
    include(joinpath("MSA", "Residues.jl"))
    include(joinpath("MSA", "Alphabet.jl"))
    include(joinpath("MSA", "ThreeLetters.jl"))
    include(joinpath("MSA", "Annotations.jl"))
    include(joinpath("MSA", "MultipleSequenceAlignment.jl"))
    include(joinpath("MSA", "GeneralParserMethods.jl"))
    include(joinpath("MSA", "IO.jl"))
    include(joinpath("MSA", "General.jl"))
    include(joinpath("MSA", "MSAEditing.jl"))
    include(joinpath("MSA", "MSAStats.jl"))
    include(joinpath("MSA", "Shuffle.jl"))
    include(joinpath("MSA", "Identity.jl"))
    include(joinpath("MSA", "Hobohm.jl"))
    include(joinpath("MSA", "MSAAnnotations.jl"))
end

# Information
@testset "Information" begin
    include(joinpath("Information", "ContingencyTables.jl"))
    include(joinpath("Information", "Counters.jl"))
    include(joinpath("Information", "InformationMeasures.jl"))
    include(joinpath("Information", "Iterations.jl"))
    include(joinpath("Information", "CorrectedMutualInformation.jl"))
    include(joinpath("Information", "Gaps.jl"))
    include(joinpath("Information", "Externals.jl"))
end

# PDB
@testset "PDB" begin
    include(joinpath("PDB", "PDB.jl"))
    include(joinpath("PDB", "Contacts.jl"))
    include(joinpath("PDB", "Kabsch.jl"))
end

# SIFTS
@testset "SIFTS" begin
    include(joinpath("SIFTS", "SIFTS.jl"))
end

# Pfam
@testset "Pfam" begin
    include(joinpath("Pfam", "Pfam.jl"))
end

# # Scripts
# include("scripts.jl")

print("""

----- =D -----

""")
