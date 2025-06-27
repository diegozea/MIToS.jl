# TODO: Identify chain-breaks
# I like this for the peptide bond: https://www.ruppweb.org/Xray/tutorial/engtabl.html
# For Calpha-Calpha: 
#  "The frequency distribution of the distance of consecutive Ca atoms
# in ~100 proteins in the top100H database (a database consisting of
# high quality structures) 27 shows that the distance between consecutive Ca atoms are 
# distributed normally with a mean of 3.8 Å and standard deviation of 0.04 Å (Figure 1a). 
# Out of 16,162 pairs of consecutive Ca atom distances, 14,281 (88%) were spaced 3.8 Å
# apart, 1297 (8%) were spaced 3.9 Å apart and 553 (3%) were spaced 3.7 Å apart. 
# Only 31 (0.1%) pairs of consecutive Ca atom distances had values different than these 
# (highest being 4 Å and the lowest being 2.9 Å)." https://f1000researchdata.s3.amazonaws.com/manuscripts/3062/f7f90715-66e7-450a-8149-87ca280b8afd_2148%20-%20sandeep%20chakraborty%20v2.pdf


# Plot coordinates of the C alpha with best occupancy
@recipe function plot(residues::AbstractVector{PDBResidue})
    ca = CAmatrix(residues)
    chains = [r.id.chain for r in residues if r.id.group == "ATOM"]
    group --> chains
    ca[:, 1], ca[:, 2], ca[:, 3]
end
