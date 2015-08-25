@auto_hash_equals immutable PDBResidueIdentifier
  PDBe_number::Nullable{Int} # PDBe
	PDB_number::ASCIIString # PDB
	name::ASCIIString
	group::ASCIIString
	model::ASCIIString
	chain::ASCIIString
end

# function !=(a::PDBResidueIdentifier, b::PDBResidueIdentifier)
#   a.number != b.number && a.name != b.name && a.group != b.group && a.chain != b.chain && a.model != b.model
# end

@auto_hash_equals immutable Coordinates{T<:AbstractFloat}
  x::T
  y::T
  z::T
end

distance(a::Coordinates, b::Coordinates) = sqrt((a.x - b.x)^2 + (a.y - b.y)^2 + (a.z - b.z)^2)

contact(a::Coordinates, b::Coordinates, limit::AbstractFloat) = distance(a,b) <= limit ? true : false

-(a::Coordinates, b::Coordinates) = Coordinates(a.x - b.x, a.y - b.y, a.z - b.z)

norm(a::Coordinates) = sqrt(a.x^2 + a.y^2 + a.z^2)

dot(a::Coordinates, b::Coordinates) = (a.x * b.x) + (a.y * b.y) + (a.z * b.z)

"""Angle (in degrees) at b between a-b and b-c"""
function angle(a::Coordinates, b::Coordinates, c::Coordinates)
  A = b - a
  B = b - c
  norms = (norm(A)*norm(B))
  if norms != 0
    return( acosd(dot(A,B) / norms) )
  else
    return(0.0)
  end
end

vec{T}(a::Coordinates{T}) = T[a.x, a.y, a.z]

function cross(a::Coordinates, b::Coordinates)
  normal = cross(vec(a), vec(b))
  Coordinates(normal[1], normal[2], normal[3])
end

@auto_hash_equals immutable PDBAtom{T<:AbstractFloat}
  residueid::PDBResidueIdentifier
  coordinates::Coordinates{T}
  atomid::ASCIIString
  element::ASCIIString
  occupancy::T
  B::ASCIIString
end

distance(a::PDBAtom, b::PDBAtom) = distance(a.coordinates, b.coordinates)

contact(a::PDBAtom, b::PDBAtom, limit::AbstractFloat) = contact(a.coordinates, b.coordinates, limit)

angle(a::PDBAtom, b::PDBAtom, c::PDBAtom) = angle(a.coordinates, b.coordinates, c.coordinates)

cross(a::PDBAtom, b::PDBAtom) = cross(a.coordinates, b.coordinates)

type PDBResidue{T<:AbstractFloat}
  id::PDBResidueIdentifier
	atoms::Vector{PDBAtom{T}}
end

length(res::PDBResidue) = length(res.atoms)

function findheavy(res::PDBResidue)
  N = length(res)
  indices = Array(Int,N)
  j = 0
  @inbounds for i in 1:N
    if res.atoms[i].element != "H"
      j += 1
      indices[j] = i
    end
  end
  resize!(indices, j)
end

function findatoms(res::PDBResidue, atomid::ASCIIString)
  N = length(res)
  indices = Array(Int,N)
  j = 0
  @inbounds for i in 1:N
    if res.atoms[i].atomid == atomid
      j += 1
      indices[j] = i
    end
  end
  resize!(indices, j)
end

function findCB(res::PDBResidue)
  N = length(res)
  indices = Array(Int,N)
  j = 0
  @inbounds for i in 1:N
    if (res.atoms[i].residueid.name == "GLY" && res.atoms[i].atomid == "CA") || res.atoms[i].atomid == "CB"
      j += 1
      indices[j] = i
    end
  end
  resize!(indices, j)
end

function selectbestoccupancy(res::PDBResidue, indices::Vector{Int})
  Ni = length(indices)
  if Ni == 1
    return(indices[1])
  end
  Na = length(res)
  if Ni == 0 || Ni > Na
    throw("There are not atom indices or they are more atom indices than atoms in the Residue")
  end
  indice = 1
  occupancy = 0.0
  for i in indices
    actual_occupancy = res.atoms[i].occupancy
    if actual_occupancy > occupancy
      occupancy = actual_occupancy
      indice = i
    end
  end
  return(indice)
end

function _update_distance(a, b, i, j, dist)
  actual_dist = distance(a.atoms[i], b.atoms[j])
  if actual_dist < dist
    return(actual_dist)
  else
    return(dist)
  end
end

"""Heavy, All, CA, CB (CA for GLY)"""
function distance(a::PDBResidue, b::PDBResidue; criteria::ASCIIString="All")
  dist = Inf
  if criteria == "All"
    Na = length(a)
    Nb = length(b)
    @inbounds for i in 1:Na
      for j in 1:Nb
        dist = _update_distance(a, b, i, j, dist)
      end
    end
  elseif criteria == "Heavy"
    indices_a = findheavy(a)
    indices_b = findheavy(b)
    if length(indices_a) != 0 && length(indices_b) != 0
      for i in indices_a
        for j in indices_b
          dist = _update_distance(a, b, i, j, dist)
        end
      end
    end
  elseif criteria == "CA"
    indices_a = findatoms(a, "CA")
    indices_b = findatoms(b, "CA")
    if length(indices_a) != 0 && length(indices_b) != 0
      for i in indices_a
        for j in indices_b
          dist = _update_distance(a, b, i, j, dist)
        end
      end
    end
  elseif criteria == "CB"
    indices_a = findCB(a)
    indices_b = findCB(b)
    if length(indices_a) != 0 && length(indices_b) != 0
      for i in indices_a
        for j in indices_b
          dist = _update_distance(a, b, i, j, dist)
        end
      end
    end
  end
  dist
end

"""Test if the function f is true for any pair of atoms between the residues a and b"""
function any(f::Function, a::PDBResidue, b::PDBResidue)
  Na = length(a)
  Nb = length(b)
  @inbounds for i in 1:Na
    for j in 1:Nb
      if f(a.atoms[i], b.atoms[j])
        return(true)
      end
    end
  end
  return(false)
end

"""Heavy, All, CA, CB (CA for GLY)"""
function contact(a::PDBResidue, b::PDBResidue, limit::AbstractFloat; criteria::ASCIIString="All")
  if criteria == "All"
    Na = length(a)
    Nb = length(b)
    @inbounds for i in 1:Na
      for j in 1:Nb
        if contact(a.atoms[i], b.atoms[j], limit)
          return(true)
        end
      end
    end
  elseif criteria == "Heavy"
    indices_a = findheavy(a)
    indices_b = findheavy(b)
    if length(indices_a) != 0 && length(indices_b) != 0
      for i in indices_a
        for j in indices_b
          if contact(a.atoms[i], b.atoms[j], limit)
            return(true)
          end
        end
      end
    end
  elseif criteria == "CA"
    indices_a = findatoms(a, "CA")
    indices_b = findatoms(b, "CA")
    if length(indices_a) != 0 && length(indices_b) != 0
      for i in indices_a
        for j in indices_b
          if contact(a.atoms[i], b.atoms[j], limit)
            return(true)
          end
        end
      end
    end
  elseif criteria == "CB"
    indices_a = findCB(a)
    indices_b = findCB(b)
    if length(indices_a) != 0 && length(indices_b) != 0
      for i in indices_a
        for j in indices_b
          if contact(a.atoms[i], b.atoms[j], limit)
            return(true)
          end
        end
      end
    end
  end
  false
end
