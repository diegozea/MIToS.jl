using Documenter

using MIToS
using MIToS.Utils
using MIToS.MSA
using MIToS.Information
using MIToS.PDB
using MIToS.SIFTS
using MIToS.Pfam
using Aqua
using LinearAlgebra
using Random
using OrderedCollections    # OrderedDict
using Statistics            # mean
using DelimitedFiles        # readdlm
using ROCAnalysis           # AUC
using Clustering            # test/MSA/Hobohm.jl
using NamedArrays           # array
using StatsBase             # WeightVec
using PairwiseListMatrices  # getlist

const DATA = joinpath(@__DIR__, "data")

@testset verbose = true "Aqua" begin
    # The ambiguities are not caused by MIToS
    Aqua.test_all(MIToS, ambiguities = false)
end

# Utils
@testset verbose = true "Utils" begin
    include("Utils/GeneralUtils.jl")
end

# MSA
@testset verbose = true "MSA" begin
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
    include("MSA/Sequences.jl")
    include("MSA/Shuffle.jl")
    include("MSA/Identity.jl")
    include("MSA/Hobohm.jl")
    include("MSA/MSAAnnotations.jl")
    include("MSA/GetIndex.jl")
    include("MSA/Concatenation.jl")
end

# Information
@testset verbose = true "Information" begin
    include("Information/ContingencyTables.jl")
    include("Information/Counters.jl")
    include("Information/InformationMeasures.jl")
    include("Information/Iterations.jl")
    include("Information/CorrectedMutualInformation.jl")
    include("Information/Gaps.jl")
    include("Information/Externals.jl")
end

# PDB
@testset verbose = true "PDB" begin
    include("PDB/PDB.jl")
    include("PDB/Contacts.jl")
    include("PDB/Kabsch.jl")
    include("PDB/Internals.jl")
    include("PDB/Sequences.jl")
    include("PDB/AlphaFoldDB.jl")
end

# SIFTS
@testset verbose = true "SIFTS" begin
    include("SIFTS/SIFTS.jl")
end

# Pfam
@testset verbose = true "Pfam" begin
    include("Pfam/Pfam.jl")
end
