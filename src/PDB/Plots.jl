# Plot coordinates of the C alpha with best occupancy
@recipe function plot(residues::AbstractVector{PDBResidue})
    ca = CAmatrix(residues)
    ca[:,1], ca[:,2], ca[:,3]
end
