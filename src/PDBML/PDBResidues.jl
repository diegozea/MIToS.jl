type PDBResidue
	number::ASCIIString
	name::ASCIIString
	group::ASCIIString
	atom::Vector{ASCIIString}
	coordinates::Vector{(FloatingPoint,FloatingPoint,FloatingPoint)}
	element::Vector{ASCIIString}
	occupancy::Vector{FloatingPoint}
	B::Vector{ASCIIString}
end

function pdbresidue(number::ASCIIString, name::ASCIIString,
	group::ASCIIString, atom::ASCIIString,
	coordinates::(FloatingPoint,FloatingPoint,FloatingPoint), element::ASCIIString,
	occupancy::FloatingPoint, B::ASCIIString)
	PDBResidue(number, name, group, [atom], [coordinates], [element], [occupancy], [B])
end

function __add_atom!(res::PDBResidue, number::ASCIIString, name::ASCIIString,
	group::ASCIIString, atom::ASCIIString,
	coordinates::(FloatingPoint,FloatingPoint,FloatingPoint), element::ASCIIString,
	occupancy::FloatingPoint, B::ASCIIString)
	if res.number == number && res.name == name && res.group == group
		push!(res.coordinates, coordinates)
		push!(res.atom, atom)
		push!(res.element, element)
		push!(res.occupancy, occupancy)
		push!(res.B, B)
		return(res)
	else
		throw("It isn't the same residue: $((res.number, res.name, res.group)) != $((number, name, group)) ")
	end
end
