using MIToS
using Base.Test

# Utils
include("utils.jl")
include("eachlinestring.jl")
# MSA
include("residues.jl")
include("indexedarrays.jl")
include("annotations.jl")
include("rawalnandgaps.jl")
include("multiplesequencealignment.jl")
include("msaannotations.jl")
include("shuffle.jl")
# Clustering
include("clustering.jl")
# PDB
include("pdb.jl")
include("contacts.jl")
#  SIFTS
include("sifts.jl")
# Information
include("probabilities.jl")
include("informationmeasures.jl")
include("iterations.jl")
include("buslje09.jl")
