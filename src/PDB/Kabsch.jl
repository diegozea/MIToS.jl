# Jorge Fernández de Cossío Díaz ( @cossio ) wrote the kabsch, rmsd and center! functions.

"""
`kabsch(A::Matrix{Float64}, B::Matrix{Float64})`

This function takes two sets of points, `A` (refrence) and `B` as NxD matrices, where D is the dimension and N is the number of points.
Assumes that the centroids of `A` and `B` are at the origin of coordinates.
You can call `center!` on each matrix before calling `kabsch` to center the matrices in the `(0.0, 0.0, 0.0)`.
Rotates `B` so that `rmsd(A,B)` is minimized.
Returns the rotation matrix. You should do `B * RotationMatrix` to get the rotated B.
"""
function kabsch(A::Matrix{Float64}, B::Matrix{Float64})
    @assert size(A) == size(B)
    M::Matrix{Float64} = B' * A
    χ = eye(M)
    χ[end,end] = sign(det(M))
    u::Matrix{Float64}, σ::Vector{Float64}, v::Matrix{Float64} = svd(M)
    return u * χ * v'
end

"""
`center!(A::Matrix{Float64})`

Takes a set of points `A` as an NxD matrix (N: number of points, D: dimension).
Translates `A` in place so that its centroid is at the origin of coordinates
"""
function center!(A::Matrix{Float64})
    for i = 1:size(A)[2]
        A[:,i] -= mean(A[:,i])
    end
end

"""
`rmsd(A::Matrix{Float64}, B::Matrix{Float64})`

Return RMSD between two sets of points `A` and `B`, given as NxD matrices (N: number of points, D: dimension).
"""
function rmsd(A::Matrix{Float64}, B::Matrix{Float64})
    @assert size(A) == size(B)
    N, D = size(A)
    s::Float64 = 0.0
    for i=1:N, j=1:D
        s += (B[i,j]::Float64 - A[i,j]::Float64)^2
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
function coordinatesmatrix(res::PDBResidue)
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

coordinatesmatrix(residues::AbstractVector{PDBResidue}) = vcat(map(coordinatesmatrix, residues)...)

"""
Returns a `Matrix{Float64}` with the centered coordinates of all the atoms in `residues`.
An optional positional argument `CA` (default: `true`) defines if only Cα carbons should be used to center the matrix.
"""
function centeredcoordinates(residues::AbstractVector{PDBResidue}, CA::Bool=true)
    coordinates = PDB.coordinatesmatrix(residues)
    meancoord = CA ? mean(PDB.CAmatrix(residues), 1) : mean(coordinates, 1)
    coordinates .- meancoord
end

"""
Returns a new `Vector{PDBResidue}` with the `PDBResidue`s having centered coordinates.
An optional positional argument `CA` (default: `true`) defines if only Cα carbons should be used to center the matrix.
"""
function centeredresidues(residues::AbstractVector{PDBResidue}, CA::Bool=true)
    coordinates = centeredcoordinates(residues, CA)
    _change_coordinates(residues, coordinates)
end

"Returns a new `PDBAtom` but with a new `coordinates`"
function _change_coordinates(atom::PDBAtom, coordinates::Coordinates)
    PDBAtom(coordinates,
            copy(atom.atom),
            copy(atom.element),
            copy(atom.occupancy),
            copy(atom.B))
end

"Returns a new `Vector{PDBResidues}` with (x,y,z) from a coordinates `Matrix{Float64}`"
function _change_coordinates(residues::AbstractVector{PDBResidue}, coordinates::Matrix{Float64})
    nres = length(residues)
    updated = Array(PDBResidue, nres)
    j = 0
    for i in 1:nres
        res = residues[i]
        centeredatoms = map(res.atoms) do atom
            _change_coordinates(atom, Coordinates(vec(coordinates[j += 1,:])))
        end
        updated[i] = PDBResidue(res.id, centeredatoms)
    end
    updated
end

"Returns a new `PDBAtom` but with a `B` as B-factor"
function _change_B(atom::PDBAtom, B::ASCIIString)
    PDBAtom(copy(atom.coordinates),
            copy(atom.atom),
            copy(atom.element),
            copy(atom.occupancy),
            B)
end

_iscentered(x::Float64, y::Float64, z::Float64) = (abs(x) <= 1e-14) && (abs(y) <= 1e-14) && (abs(z) <= 1e-14)

_iscentered(meanCα::Vector{Float64}) = _iscentered(meanCα[1], meanCα[2], meanCα[3])

_iscentered(CA::Matrix{Float64}) = _iscentered(vec(mean(CA,1)))

"""
This function takes `A::Vector{PDBResidue}` (reference) and `B::Vector{PDBResidue}`.
Translates `A` and `B` to the origin of coordinates,
and rotates `B` so that `rmsd(A,B)` is minimized with the Kabsch algorithm (using only their α carbons).
Returns the rotated and translated versions of `A` and `B`, and the RMSD value.
"""
function superimpose(A::Vector{PDBResidue}, B::Vector{PDBResidue})
    @assert length(A) == length(B)
    Bxyz = PDB.coordinatesmatrix(B)
    ACα = PDB.CAmatrix(A)
    BCα = PDB.CAmatrix(B)
    meanACα = vec(mean(ACα,1))
    meanBCα = vec(mean(BCα,1))
    if !_iscentered(meanBCα)
        @inbounds BCα[:,1] -= meanBCα[1]
        @inbounds BCα[:,2] -= meanBCα[2]
        @inbounds BCα[:,3] -= meanBCα[3]
        @inbounds Bxyz[:,1] -= meanBCα[1]
        @inbounds Bxyz[:,2] -= meanBCα[2]
        @inbounds Bxyz[:,3] -= meanBCα[3]
    end
    if !_iscentered(meanACα)
        @inbounds ACα[:,1] -= meanACα[1]
        @inbounds ACα[:,2] -= meanACα[2]
        @inbounds ACα[:,3] -= meanACα[3]
        Axyz = PDB.coordinatesmatrix(A)
        @inbounds Axyz[:,1] -= meanACα[1]
        @inbounds Axyz[:,2] -= meanACα[2]
        @inbounds Axyz[:,3] -= meanACα[3]
        RotationMatrix = PDB.kabsch(ACα, BCα)
        return(
            _change_coordinates(A, Axyz),
            _change_coordinates(B, Bxyz * RotationMatrix),
            PDB.rmsd(ACα, BCα * RotationMatrix)
            )
    else
        RotationMatrix = PDB.kabsch(ACα, BCα)
        return(
            A,
            _change_coordinates(B, Bxyz * RotationMatrix),
            PDB.rmsd(ACα, BCα * RotationMatrix)
            )
    end
end
