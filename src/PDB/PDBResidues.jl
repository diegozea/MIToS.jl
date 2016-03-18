# PDB Types
# =========

@auto_hash_equals immutable PDBResidueIdentifier
  PDBe_number::ASCIIString # PDBe
  number::ASCIIString # PDB
  name::ASCIIString
  group::ASCIIString
  model::ASCIIString
  chain::ASCIIString
end

@auto_hash_equals immutable Coordinates # <: FixedVectorNoTuple{Float64, 3}
  x::Float64
  y::Float64
  z::Float64
end

@auto_hash_equals immutable PDBAtom
  coordinates::Coordinates
  atom::ASCIIString
  element::ASCIIString
  occupancy::Float64
  B::ASCIIString
end

type PDBResidue
  id::PDBResidueIdentifier
  atoms::Vector{PDBAtom}
end

length(res::PDBResidue) = length(res.atoms)

# Distances and geom
# ==================

distance(a::Coordinates, b::Coordinates) = sqrt((a.x - b.x)^2 + (a.y - b.y)^2 + (a.z - b.z)^2)

"distance(a,b) <= limit"
contact(a::Coordinates, b::Coordinates, limit::AbstractFloat) = distance(a,b) <= limit

#-------------------------------------------------------------------------------#
vec(a::Coordinates) = Float64[a.x, a.y, a.z]
norm(a::Coordinates) = sqrt(a.x^2 + a.y^2 + a.z^2)
dot(a::Coordinates, b::Coordinates) = (a.x * b.x) + (a.y * b.y) + (a.z * b.z)
size(a::Coordinates) = (3,)
length(a::Coordinates) = 3
-(a::Coordinates, b::Coordinates) = Coordinates(a.x - b.x, a.y - b.y, a.z - b.z)
+(a::Coordinates, b::Coordinates) = Coordinates(a.x + b.x, a.y + b.y, a.z + b.z)
./(a::Coordinates, b::Int) = Coordinates(a.x / b, a.y / b, a.z / b )
function cross(a::Coordinates, b::Coordinates)
  normal = cross(vec(a), vec(b))
  Coordinates(normal[1], normal[2], normal[3])
end
#-------------------------------------------------------------------------------#

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

distance(a::PDBAtom, b::PDBAtom) = distance(a.coordinates, b.coordinates)

contact(a::PDBAtom, b::PDBAtom, limit::Float64) = contact(a.coordinates, b.coordinates, limit)

angle(a::PDBAtom, b::PDBAtom, c::PDBAtom) = angle(a.coordinates, b.coordinates, c.coordinates)

cross(a::PDBAtom, b::PDBAtom) = cross(a.coordinates, b.coordinates)

# Find Residues/Atoms
# ===================

isobject(res::PDBResidue, tests::AbstractTest...) = isobject(res.id, tests...)

findobjects(res::PDBResidue, tests::AbstractTest...) = findobjects(res.atoms, tests...)

# @residues
# =========

_test_stringfield(field::Symbol, test::Real) = Is(field, string(test))
_test_stringfield(field::Symbol, test::Union{Char, Symbol}) = Is(field, string(test))
_test_stringfield(field::Symbol, test::Union{ASCIIString, Regex, Function}) = Is(field, test)
_test_stringfield(field::Symbol, test::Union{UnitRange, IntSet, Set, Array, Base.KeyIterator}) = In(field, test)

_is_wildcard(test::ASCIIString) = test == "*"
_is_wildcard(test::Char) = test == '*'
_is_wildcard(test) = false

function _add_test_stringfield!(tests::Vector{TestType}, field::Symbol, test)
  if !_is_wildcard(test)
    push!(tests, _test_stringfield(field, test))
  end
  tests
end

function _residues_tests(model, chain, group, residue)
  args = TestType[]
  _add_test_stringfield!(args, :model, model)
  _add_test_stringfield!(args, :chain, chain)
  _add_test_stringfield!(args, :group, group)
  _add_test_stringfield!(args, :number, residue)
  args
end

residues(residue_list, model, chain, group, residue) = collectobjects(residue_list, _residues_tests(model, chain, group, residue)...)

macro residues(residue_list,
               model::Symbol, m, #::Union(Int, Char, ASCIIString, Symbol),
               chain::Symbol, c, #::Union(Char, ASCIIString, Symbol),
               group::Symbol, g,
               residue::Symbol, r)
  if model == :model && chain == :chain && group == :group && residue == :residue
    return :(residues($(esc(residue_list)), $(esc(m)), $(esc(c)), $(esc(g)), $(esc(r))))
  else
    throw(ArgumentError("The signature is @residues ___ model ___ chain ___ group ___ residue ___"))
  end
end

function residuesdict(residue_list, model, chain, group, residue)
  res_index = findobjects(residue_list, _residues_tests(model, chain, group, residue)...)
  dict = sizehint!( OrderedDict{ASCIIString, PDBResidue}() , length(res_index) )
  for i in res_index
    dict[ residue_list[i].id.number ] = residue_list[i]
  end
  dict
end

macro residuesdict(residue_list,
                   model::Symbol, m, #::Union(Int, Char, ASCIIString, Symbol),
                   chain::Symbol, c, #::Union(Char, ASCIIString, Symbol),
                   group::Symbol, g,
                   residue::Symbol, r)
  if model == :model && chain == :chain && group == :group && residue == :residue
    return :(residuesdict($(esc(residue_list)), $(esc(m)), $(esc(c)), $(esc(g)), $(esc(r))))
  else
    throw(ArgumentError("The signature is @residues ___ model ___ chain ___ group ___ residue ___"))
  end
end

# @atoms
# ======

function atoms(residue_list, model, chain, group, residue, atom)
  _is_wildcard(atom) ?
    vcat(Vector{PDBAtom}[ res.atoms for res in residues(residue_list, model, chain, group, residue) ]...) :
    vcat(Vector{PDBAtom}[ collectobjects(res.atoms, _test_stringfield(:atom, atom)) for res in residues(residue_list, model, chain, group, residue) ]...)
end

macro atoms(residue_list,
               model::Symbol, m, #::Union(Int, Char, ASCIIString, Symbol),
               chain::Symbol, c, #::Union(Char, ASCIIString, Symbol),
               group::Symbol, g,
               residue::Symbol, r,
               atom::Symbol, a)
  if model == :model && chain == :chain && group == :group && residue == :residue && atom == :atom
    return :(atoms($(esc(residue_list)), $(esc(m)), $(esc(c)), $(esc(g)), $(esc(r)), $(esc(a))))
  else
    throw(ArgumentError("The signature is @atoms ___ model ___ chain ___ group ___ residue ___ atom ___"))
  end
end

# Special find...
# ===============

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

function findatoms(res::PDBResidue, atom::ASCIIString)
  N = length(res)
  indices = Array(Int,N)
  j = 0
  @inbounds for i in 1:N
    if res.atoms[i].atom == atom
      j += 1
      indices[j] = i
    end
  end
  resize!(indices, j)
end

"Returns a vector of indices for CB (CA for GLY)"
function findCB(res::PDBResidue)
  N = length(res)
  indices = Array(Int,N)
  atom = res.id.name == "GLY" ? "CA" : "CB"
  j = 0
  @inbounds for i in 1:N
    if res.atoms[i].atom == atom
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

# More Distances
# ==============

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

# For Aromatic
# ============

function _get_plane(residue::PDBResidue)
  name = residue.id.name
  planes = Vector{PDBAtom}[]
  if name != "TRP"
    plane = PDBAtom[]
    for atom in residue.atoms
      if (name, atom.atom) in PDB._aromatic
        push!(plane, atom)
      end
    end
    push!(planes, plane)
  else
    plane1 = PDBAtom[]
    plane2 = PDBAtom[]
    for atom in residue.atoms
      if atom.atom in Set{ASCIIString}(["CE2", "CD2", "CZ2", "CZ3", "CH2", "CE3"])
        push!(plane1, atom)
      end
      if atom.atom in Set{ASCIIString}(["CG", "CD1", "NE1", "CE2", "CD2"])
        push!(plane2, atom)
      end
    end
    push!(planes, plane1)
    push!(planes, plane2)
  end
  planes
end

function bestoccupancy!(atoms::Vector{PDBAtom})
  N = length(atoms)
  name = atoms[1].atom
  occupancy = atoms[1].occupancy
  for i in 2:N
    if atoms[i].occupancy > occupancy && atoms[i].atom == name
      splice!(atoms, i-1)
    end
    name = atoms[i].atom
    occupancy = atoms[i].occupancy
  end
  atoms
end

function _centre(planes::Vector{Vector{PDBAtom}})
  subset = Vector{PDBAtom}[ bestoccupancy!(atoms) for atoms in planes ]
  polyg = Vector{Coordinates}[ Coordinates[ a.coordinates for a in atoms ] for atoms in subset ]
  Coordinates[ sum(points)./length(points) for  points in polyg ]
end

# function _simple_normal_and_centre(atoms::Vector{PDBAtom})
#   bestoccupancy!(atoms)
#   points = Coordinates[ a.coordinates for a in atoms ]
#   (cross(points[2] - points[1], points[3] - points[1]), sum(points)./length(points))
# end

# For Parsing...
# ==============

"Used for parsing a PDB file into Vector{PDBResidue}"
function _generate_residues(residue_dict::OrderedDict{PDBResidueIdentifier, Vector{PDBAtom}})
  len = length(residue_dict)
  residue_vect = Array(PDBResidue, len)
  i = 1
  for (k,v) in residue_dict
    residue_vect[i] = PDBResidue(k, v)
    i += 1
  end
  residue_vect
end

# Show PDB* objects (using Formatting)
# ====================================

const _Format_ResidueID = FormatExpr("{:>15} {:>15} {:>15} {:>15} {:>15} {:>15}\n")
const _Format_ATOM = FormatExpr("{:>50} {:>15} {:>15} {:>15} {:>15}\n")

const _Format_ResidueID_Res = FormatExpr("\t\t{:>15} {:>15} {:>15} {:>15} {:>15} {:>15}\n")
const _Format_ATOM_Res = FormatExpr("\t\t{:<10} {:>50} {:>15} {:>15} {:>15} {:>15}\n")

function show(io::IO, id::PDBResidueIdentifier)
  printfmt(io, _Format_ResidueID, "PDBe_number", "number", "name", "group", "model", "chain")
  printfmt(io, _Format_ResidueID, string('"',id.PDBe_number,'"'), string('"',id.number,'"'), string('"',id.name,'"'),
           string('"',id.group,'"'), string('"',id.model,'"'), string('"',id.chain,'"'))
end

function show(io::IO, atom::PDBAtom)
  printfmt(io, _Format_ATOM, "coordinates", "atom", "element", "occupancy", "B")
  printfmt(io, _Format_ATOM, atom.coordinates, string('"',atom.atom,'"'), string('"',atom.element,'"'),
           atom.occupancy, string('"',atom.B,'"'))
end

function show(io::IO, res::PDBResidue)
  println(io, "PDBResidue:\n\tid::PDBResidueIdentifier")
  printfmt(io, _Format_ResidueID_Res, "PDBe_number", "number", "name", "group", "model", "chain")
  printfmt(io, _Format_ResidueID_Res, string('"',res.id.PDBe_number,'"'), string('"',res.id.number,'"'),
           string('"',res.id.name,'"'), string('"',res.id.group,'"'), string('"',res.id.model,'"'), string('"',res.id.chain,'"'))
  len = length(res)
  println(io, "\tatoms::Vector{PDBAtom}\tlength: ", len)
  for i in 1:len
    printfmt(io, _Format_ATOM_Res, "", "coordinates", "atom", "element", "occupancy", "B")
    printfmt(io, _Format_ATOM_Res, string(i,":"), res.atoms[i].coordinates, string('"',res.atoms[i].atom,'"'),
             string('"',res.atoms[i].element,'"'), res.atoms[i].occupancy, string('"',res.atoms[i].B,'"'))
  end
end
