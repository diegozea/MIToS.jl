using Base.Test

using MIToS: Utils, MSA

using NamedArrays
using PairwiseListMatrices  # getlist

# using MIToS.Information
# using MIToS.PDB
# using MIToS.SIFTS
# using MIToS.Pfam
#

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
    include(joinpath("MSA", "ContingencyTables.jl"))
end

# # MSA
# include("indexedarrays.jl") # TO DO: Agregar si es utilizado luego
# # PDB
# include("pdb.jl")
# include("contacts.jl")
# include("kabsch.jl")
# #  SIFTS
# include("sifts.jl")
# # Information
# include("probabilities.jl")
# include("informationmeasures.jl")
# include("iterations.jl")
# include("buslje09.jl")
# # Pfam
# include("pfam.jl")
# # Scripts
# include("scripts.jl")

print("""

----- =D -----

""")
