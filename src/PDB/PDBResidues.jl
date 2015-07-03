import Base: ==, !=, hash, isequal

immutable PDBResidueIdentifier
	number::ASCIIString
	name::ASCIIString
	group::ASCIIString
	model::ASCIIString
	chain::ASCIIString
end

hash(a::PDBResidueIdentifier) = hash(string(a.number, a.name, a.group, a.model,a.chain))

isequal(a::PDBResidueIdentifier, b::PDBResidueIdentifier) = hash(a) == hash(b)

function ==(a::PDBResidueIdentifier, b::PDBResidueIdentifier)
  a.number == b.number && a.name == b.name && a.group == b.group && a.chain == b.chain && a.model == b.model
end

function !=(a::PDBResidueIdentifier, b::PDBResidueIdentifier)
  a.number != b.number && a.name != b.name && a.group != b.group && a.chain != b.chain && a.model != b.model
end

immutable Coordinates{T<:FloatingPoint}
  x::T
  y::T
  z::T
end

distance(a::Coordinates, b::Coordinates) = sqrt((a.x - b.x)^2 + (a.y - b.y)^2 + (a.z - b.z)^2)

contact(a::Coordinates, b::Coordinates, limit::FloatingPoint) = distance(a,b) <= limit ? true : false

immutable PDBAtom{T<:FloatingPoint}
  residueid::PDBResidueIdentifier
  coordinates::Coordinates{T}
  atomid::ASCIIString
  element::ASCIIString
  occupancy::T
  B::ASCIIString
end

distance(a::PDBAtom, b::PDBAtom) = distance(a.coordinates, b.coordinates)

contact(a::PDBAtom, b::PDBAtom, limit::FloatingPoint) = contact(a.coordinates, b.coordinates, limit)

type PDBResidue{T<:FloatingPoint}
  id::PDBResidueIdentifier
	atoms::Vector{PDBAtom{T}}
end
