using Base.Test
using MIToS
using MIToS.Utils
using MIToS.MSA
# using MIToS.Information
# using MIToS.PDB
# using MIToS.SIFTS
# using MIToS.Pfam
# using PairwiseListMatrices

const DATA = joinpath(pwd(), "data")

# Utils
@testset "Utils" begin
    include(joinpath("Utils", "GeneralUtils.jl"))
end

# MSA
@testset "MSA" begin
    include(joinpath("MSA", "Residues.jl"))
    include(joinpath("MSA", "ThreeLetters.jl"))
end

# # MSA
# include("residues.jl")
# include("indexedarrays.jl")
# include("annotations.jl")
# include("rawalnandgaps.jl")
# include("multiplesequencealignment.jl")
# include("msaannotations.jl")
# include("shuffle.jl")
# # Clustering
# include("clustering.jl")
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
