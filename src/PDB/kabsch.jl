# Jorge Fernández de Cossío Díaz ( @cossio ) wrote the kabsch, rmsd and center! functions.

"""
`kabsch(X::Matrix{Float64}, Y::Matrix{Float64})`

Input: Two sets of points X,Y as NxD matrices, where D is dimension, N number of points, and a preallocated NxD matrix rotatedX.
Assumes that the centroids of X and Y are at the origin of coordinates. Rotates X so that RMSD(X,Y) is minimized.
Output: Returns the rotation matrix. You should do `X * RotationMatrix` to get the rotated X.
"""
function kabsch(X::Matrix{Float64}, Y::Matrix{Float64})
    @assert size(X) == size(Y)
    A::Matrix{Float64} = X' * Y
    χ = eye(A)
    χ[end,end] = sign(det(A))
    u::Matrix{Float64}, σ::Vector{Float64}, v::Matrix{Float64} = svd(A)
    return u * χ * v'
end
"""
`center!(X::Matrix{Float64})`

Input: A set of points X as an NxD matrix (N: number of points, D: dimension).
Translates X in place so that its centroid is at the origin of coordinates
"""
function center!(X::Matrix{Float64})
    for i = 1:size(X)[2]
        X[:,i] -= mean(X[:,i])
    end
end

"""
`rmsd(A::Matrix{Float64}, B::Matrix{Float64})`

Return RMSD between two sets of points A and B, given as NxD matrices (N: number of points, D: dimension).
"""
function rmsd(A::Matrix{Float64}, B::Matrix{Float64})
    @assert size(A) == size(B)
    N, D = size(A)
    s::Float64 = 0.0
    for i=1:N, j=1:D
        s += (A[i,j]::Float64 - B[i,j]::Float64)^2
    end
    return sqrt(s / N)
end

# Kabsch for Vector{PDBResidues}
# ==============================

"Returns the Cα with best occupancy in the `PDBResidue`."
function getCA(res::PDBResidue)
    @assert length(res) != 0
    caindex = selectbestoccupancy(res, findatoms(res, "CA"))
    res.atoms[caindex]
end

"Returns a matrix with the x, y and z coordinates of the Cα with best occupancy for each `PDBResidue`."
function CAmatrix(residues::AbstractVector{PDBResidue})
    len = length(residues)
    CAmat = Array(Float64, 3, len)
    for i in 1:len
        coord = getCA(residues[i]).coordinates
        CAmat[1,i] = coord.x
        CAmat[2,i] = coord.y
        CAmat[3,i] = coord.z
    end
    CAmat'
end

"Returns a matrix with the x, y, z coordinates of each atom in each `PDBResidue`"
function coordinatematrix(res::PDBResidue)
    atoms = res.atoms
    len = length(atoms)
    mat = Array(Float64, 3, len)
    for i in 1:len
        coord = atoms[i].coordinates
        mat[1,i] = coord.x
        mat[2,i] = coord.y
        mat[3,i] = coord.z
    end
    mat'
end

coordinatematrix(residues::AbstractVector{PDBResidue}) = vcat(map(coordinatematrix, residues)...)



# a = [ 1  1 0
#       2  1 0
#     1+cos(pi/4) 1-sin(pi/4) 0 ]
# #      1  0 0 ]
# b = [  1 1 0
# #    1 2 0
#     1+cos(pi/4) 1+sin(pi/4)  0 ]
# scatter3d(a[:,1], a[:,2], a[:,3])
# scatter3d!(b[:,1], b[:,2], b[:,3])
# center!(b)
# sa = a[1:2,:]
# #sa = a[1:3,:]
# ma = mean(sa,1)
# sa = sa .- ma
# a  = a  .- ma
# scatter3d(sa[:,1], sa[:,2], sa[:,3])
# scatter3d!(b[:,1], b[:,2], b[:,3])
# F = kabsch(sa, b)
# R = a * F
# scatter3d(R[:,1], R[:,2], R[:,3], alpha=0.5)
# scatter3d!(b[:,1], b[:,2], b[:,3], alpha=0.5)
# scatter(R[:,1:2], alpha=0.75, legend=false)
# scatter!(b[:,1:2], alpha=0.5)


