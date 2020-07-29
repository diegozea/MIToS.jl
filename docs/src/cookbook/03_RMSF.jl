# # Linking structural and evolutionary information
#
#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__cookbook/notebooks/03_RMSF.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__cookbook/notebooks/03_RMSF.ipynb)
#
#
# ## Problem description
#
# The [Root Mean Squared Fluctuation (RMSF)](https://en.wikipedia.org/wiki/Mean_squared_displacement) 
# is a common way to measure residue flexibility in a structural ensemble. 
# It is a measure of how far is the residue moving from its average position 
# in the group of structures. Usually, we represent a residue position with 
# the spatial coordinates of its alpha carbon. 
#
# The protein structures should be previously superimposed to calculate the 
# RMSF, for example, by using the `superimpose` function of the 
# [`PDB` module of `MIToS`](@ref Module-PDB). In this example, we are going 
# to measure the RMSF of each residue from an NMR ensemble using the 
# `rmsf` function. 
#
# The structure superimposition could be the most complicated step of the 
# process, depending on the input data. In particular, it structures come 
# from different PDB structures or homologous proteins can require the use 
# of external programs, 
# as [MAMMOTH-mult](https://ub.cbm.uam.es/software/online/mamothmult.php) or 
# [MUSTANG](https://lcb.infotech.monash.edu/mustang/) among others, 
# tailored for this task. 
#
# In this case, we are going to use an NMR ensemble. Therefore, we are not 
# going to need to superimpose the structures as NMR models have the 
# same protein sequence and are, usually, well-aligned.
# 
#
# ## MIToS solution

import MIToS
using MIToS.PDB
using Plots

# Lets read the NMR ensemble:


pdb_file   = abspath(pathof(MIToS), "..", "..", "test", "data", "1AS5.pdb")
pdb_res = read(pdb_file, PDBFile)

# We can get an idea of the alpha carbon positions by plotting these residues:

scatter(pdb_res, legend=false)

# As we saw in the previous plot, the structure doesn't need to be 
# superimposed. Now, we are going to separate each model into different 
# vectors, storing each vector into a `Dict`:

models = Dict{String, Vector{PDBResidue}}()
for res in pdb_res
	push!(get!(models, res.id.model, []), res) 
end

# Then, we simply need to collect all the PDB models in the values 
# of the `Dict`, to get the vector of `PDBResidues` vectors required 
# to calculate the RMSF.

pdb_models = collect(values(superimposed_models))

# And, finally, call the `rmsf` function on the list of structures:

RMSF = rmsf(pdb_models)

# This return the vector of RMSF values for each residue, calculated using 
# the coordinates of the alpha carbons.  
# You can plot this vector to get an idea of the which are the most flexible 
# position in your structure:

plot(RMSF, legend=false, xlab="Residue", ylab="RMSF [Å]")
#
