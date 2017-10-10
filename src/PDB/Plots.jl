# Plot coordinates of the C alpha with best occupancy
@recipe function plot(residues::AbstractVector{PDBResidue})
    ca = CAmatrix(residues)
    chains = [r.id.chain for r in residues if r.id.group == "ATOM"]

    seriestype  -->  :scatter
    group --> chains

    ca[:,1], ca[:,2], ca[:,3]
end
