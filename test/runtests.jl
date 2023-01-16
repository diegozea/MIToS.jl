using Test
using Documenter

using MIToS
using MIToS.Utils
using MIToS.MSA
using MIToS.Information
using MIToS.PDB
using MIToS.SIFTS
using MIToS.Pfam
using LinearAlgebra
using Random
using Statistics            # mean
using DelimitedFiles        # readdlm
using ROCAnalysis           # AUC
using Clustering            # test/MSA/Hobohm.jl
using NamedArrays           # array
using StatsBase             # WeightVec
using PairwiseListMatrices  # getlist

const DATA = "data"

# Utils
@testset "Utils" begin
    include("Utils/GeneralUtils.jl")
end

# MSA
@testset "MSA" begin
    include("MSA/Residues.jl")
    include("MSA/Alphabet.jl")
    include("MSA/ThreeLetters.jl")
    include("MSA/Annotations.jl")
    include("MSA/MultipleSequenceAlignment.jl")
    include("MSA/GeneralParserMethods.jl")
    include("MSA/IO.jl")
    include("MSA/General.jl")
    include("MSA/MSAEditing.jl")
    include("MSA/MSAStats.jl")
    include("MSA/Shuffle.jl")
    include("MSA/Identity.jl")
    include("MSA/Hobohm.jl")
    include("MSA/MSAAnnotations.jl")
    include("MSA/GetIndex.jl")
    include("MSA/Concatenation.jl")
end

# Information
@testset "Information" begin
    include("Information/ContingencyTables.jl")
    include("Information/Counters.jl")
    include("Information/InformationMeasures.jl")
    include("Information/Iterations.jl")
    include("Information/CorrectedMutualInformation.jl")
    include("Information/Gaps.jl")
    include("Information/Externals.jl")
end

# PDB
@testset "PDB" begin
    include("PDB/PDB.jl")
    include("PDB/Contacts.jl")
    include("PDB/Kabsch.jl")
    include("PDB/Internals.jl")
end

# SIFTS
@testset "SIFTS" begin
    include("SIFTS/SIFTS.jl")
end

# Pfam
@testset "Pfam" begin
    include("Pfam/Pfam.jl")
end

# Scripts
@testset "Scripts" begin
    include("Scripts/Template.jl")
    include("Scripts/Scripts.jl")
end

# Doctests
doctest(MIToS)

print("""

----- =D -----

""")
