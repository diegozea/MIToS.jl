# PDB Types
# =========

"""
A `PDBResidueIdentifier` object contains the information needed to identity PDB residues.
It has the following fields that you can access at any moment for query purposes:

    - `PDBe_number` : It's only used when a PDBML is readed (PDBe number as a string).
    - `number` : PDB residue number, it includes insertion codes, e.g. `"34A"`.
    - `name` : Three letter residue name in PDB, e.g. `"LYS"`.
    - `group` : It can be `"ATOM"` or `"HETATM"`.
    - `model` : The model number as a string, e.g. `"1"`.
    - `chain` : The chain as a string, e.g. `"A"`.
"""
@auto_hash_equals struct PDBResidueIdentifier
    PDBe_number::String # PDBe
    number::String # PDB
    name::String
    group::String
    model::String
    chain::String
end

"""
A `Coordinates` object is a fixed size vector with the coordinates x,y,z.
"""
@auto_hash_equals struct Coordinates <: FieldVector{3,Float64}
    x::Float64
    y::Float64
    z::Float64
    function Coordinates(a::NTuple{3,Real})
        new(a[1], a[2], a[3])
    end
    Coordinates(x, y, z) = new(x, y, z)
end

"""
A `PDBAtom` object contains the information from a PDB atom, without information of the
residue. It has the following fields that you can access at any moment for query purposes:

    - `coordinates` : x,y,z coordinates, e.g. `Coordinates(109.641,73.162,42.7)`.
    - `atom` : Atom name, e.g. `"CA"`.
    - `element` : Element type of the atom, e.g. `"C"`.
    - `occupancy` : A float number with the occupancy, e.g. `1.0`.
    - `B` : B factor as a string, e.g. `"23.60"`.
    - `alt_id` : Alternative location ID, e.g. `"A"`.
    - `charge` : Charge of the atom, e.g. `"0"`.
"""
@auto_hash_equals struct PDBAtom
    coordinates::Coordinates
    atom::String
    element::String
    occupancy::Float64
    B::String
    alt_id::String
    charge::String
end

"""
A `PDBResidue` object contains all the information about a PDB residue. It has the
following fields that you can access at any moment for query purposes:

    - `id` : A `PDBResidueIdentifier` object.
    - `atoms` : A vector of `PDBAtom`s.
"""
@auto_hash_equals mutable struct PDBResidue
    id::PDBResidueIdentifier
    atoms::Vector{PDBAtom}
end

Base.length(res::PDBResidue) = length(res.atoms)

# Copy
# ====

# copy is not defined for String objects, using deepcopy instead

for f in (:copy, :deepcopy)

    @eval Base.$(f)(res::PDBResidueIdentifier) = PDBResidueIdentifier(
        deepcopy(res.PDBe_number),
        deepcopy(res.number),
        deepcopy(res.name),
        deepcopy(res.group),
        deepcopy(res.model),
        deepcopy(res.chain),
    )
    @eval Base.$(f)(res::Coordinates) = Coordinates($(f)(res.x), $(f)(res.y), $(f)(res.z))

    @eval Base.$(f)(res::PDBAtom) = PDBAtom(
        $(f)(res.coordinates),
        deepcopy(res.atom),
        deepcopy(res.element),
        $(f)(res.occupancy),
        deepcopy(res.B),
        deepcopy(res.alt_id),
        deepcopy(res.charge),
    )

    @eval Base.$(f)(res::PDBResidue) = PDBResidue($(f)(res.id), $(f)(res.atoms))

end

# Coordinates
# ===========

Base.vec(a::Coordinates) = Float64[a.x, a.y, a.z]

# Distances and geometry
# ----------------------

@inline function squared_distance(a::Coordinates, b::Coordinates)
    (a.x - b.x)^2 + (a.y - b.y)^2 + (a.z - b.z)^2
end

"""
It calculates the squared euclidean distance, i.e. it doesn't spend time in `sqrt`
"""
squared_distance(a::PDBAtom, b::PDBAtom) = squared_distance(a.coordinates, b.coordinates)

"""
It calculates the squared euclidean distance.
"""
distance(a::Coordinates, b::Coordinates) = sqrt(squared_distance(a, b))

distance(a::PDBAtom, b::PDBAtom) = distance(a.coordinates, b.coordinates)

function _squared_limit_contact(a::Coordinates, b::Coordinates, limit::AbstractFloat)
    squared_distance(a, b) <= limit
end

function _squared_limit_contact(a::PDBAtom, b::PDBAtom, limit::AbstractFloat)
    _squared_limit_contact(a.coordinates, b.coordinates, limit)
end

"""
`contact(a::Coordinates, b::Coordinates, limit::AbstractFloat)`

It returns true if the distance is less or equal to the limit.
It doesn't call `sqrt` because it does `squared_distance(a,b) <= limit^2`.
"""
function contact(a::Coordinates, b::Coordinates, limit::AbstractFloat)
    _squared_limit_contact(a, b, limit^2)
end

function contact(a::PDBAtom, b::PDBAtom, limit::Float64)
    contact(a.coordinates, b.coordinates, limit)
end

"""
`angle(a::Coordinates, b::Coordinates, c::Coordinates)`

Angle (in degrees) at `b` between `a-b` and `b-c`
"""
function Base.angle(a::Coordinates, b::Coordinates, c::Coordinates)
    A = b - a
    B = b - c
    norms = (norm(A) * norm(B))
    if norms != 0
        return (acosd(dot(A, B) / norms))
    else
        return (0.0)
    end
end

function Base.angle(a::PDBAtom, b::PDBAtom, c::PDBAtom)
    angle(a.coordinates, b.coordinates, c.coordinates)
end

LinearAlgebra.cross(a::PDBAtom, b::PDBAtom) = cross(a.coordinates, b.coordinates)

# Find Residues/Atoms
# ===================

"""
    ResidueQueryTypes = Union{String,Type{All},Regex,Function}

This type is used to indicate the type of the keyword arguments of functions that filter
residues or atoms, such as [`isresidue`](@ref), [`residues`](@ref), [`residuesdict`](@ref)
and [`atoms`](@ref).
"""
const ResidueQueryTypes = Union{String,Type{All},Regex,Function}

@inline _is(element::String, all::Type{All}) = true
@inline _is(element::String, value::String) = element == value
@inline _is(element::String, regex::Regex) = occursin(regex, element)
@inline _is(element::String, f::Function) = f(element)

# isresidue
# ---------

function isresidue(id::PDBResidueIdentifier, model, chain, group, residue)
    Base.depwarn(
        "isresidue using positional arguments is deprecated in favor of keyword arguments: isresidue(id; model, chain, group, residue)",
        :isresidue,
        force = true,
    )
    isresidue(id, model = model, chain = chain, group = group, residue = residue)
end

function isresidue(res::PDBResidue, model, chain, group, residue)
    Base.depwarn(
        "isresidue using positional arguments is deprecated in favor of keyword arguments: isresidue(res; model, chain, group, residue)",
        :isresidue,
        force = true,
    )
    isresidue(res.id, model = model, chain = chain, group = group, residue = residue)
end

function isresidue(
    id::PDBResidueIdentifier;
    model::ResidueQueryTypes = All,
    chain::ResidueQueryTypes = All,
    group::ResidueQueryTypes = All,
    residue::ResidueQueryTypes = All,
)
    _is(id.model, model) &&
        _is(id.chain, chain) &&
        _is(id.group, group) &&
        _is(id.number, residue)
end

"""
     isresidue(res; model=All, chain=All, group=All, residue=All)

This function tests if a `PDBResidue` has the indicated `model`, `chain`, `group`
and `residue` names/numbers. You can use the type `All` (default value) to avoid
filtering that level.
"""
function isresidue(
    res::PDBResidue;
    model::ResidueQueryTypes = All,
    chain::ResidueQueryTypes = All,
    group::ResidueQueryTypes = All,
    residue::ResidueQueryTypes = All,
)
    isresidue(res.id, model = model, chain = chain, group = group, residue = residue)
end

"""
It tests if the atom has the indicated atom name.
"""
isatom(atom::PDBAtom, name) = _is(atom.atom, name)

# select_residues
# ---------------

"""
    select_residues(residue_list; model=All, chain=All, group=All, residue=All)

This function returns a new vector with the selected subset of residues from a list of
residues. You can use the keyword arguments `model`, `chain`, `group` and `residue` to
select the residues. You can use the type `All` (default value) to avoid filtering at
a particular level.
"""
function select_residues(
    residue_list::AbstractArray{PDBResidue,N};
    model::ResidueQueryTypes = All,
    chain::ResidueQueryTypes = All,
    group::ResidueQueryTypes = All,
    residue::ResidueQueryTypes = All,
) where {N}
    filter(
        res ->
            isresidue(res, model = model, chain = chain, group = group, residue = residue),
        residue_list,
    )
end

"""
The `residues` function for `AbstractArray{PDBResidue,N}` is **deprecated**. Use the
`select_residues` function instead. So, `residues(residue_list, model, chain, group, residue)`
becomes `select_residues(residue_list; model=model, chain=chain, group=group, residue=residue)`.
"""
function residues(
    residue_list::AbstractArray{PDBResidue,N},
    model,
    chain,
    group,
    residue,
) where {N}
    Base.depwarn(
        "residues is deprecated in favor of select_residues(residue_list; model, chain, group, residue)",
        :residues,
        force = true,
    )
    select_residues(
        residue_list;
        model = model,
        chain = chain,
        group = group,
        residue = residue,
    )
end

"""
`@residues ... model ... chain ... group ... residue ...`

These return a new vector with the selected subset of residues from a list of residues. You
can use the type `All` to avoid filtering that option.

**DEPRECATED:** This macro is deprecated. Use the [`select_residues`](@ref) function instead.
"""
macro residues(
    residue_list,
    model::Symbol,
    m,
    chain::Symbol,
    c,
    group::Symbol,
    g,
    residue::Symbol,
    r,
)
    if model == :model && chain == :chain && group == :group && residue == :residue
        Base.depwarn(
            "Using the @residues macro is deprecated in favor of the select_residues function: select_residues(residue_list; model, chain, group, residue)",
            Symbol("@residues"),
            force = true,
        )
        return :(select_residues(
            $(esc(residue_list));
            model = $(esc(m)),
            chain = $(esc(c)),
            group = $(esc(g)),
            residue = $(esc(r)),
        ))
    else
        throw(
            ArgumentError(
                "The signature is @residues ___ model ___ chain ___ group ___ residue ___",
            ),
        )
    end
end

# residuesdict
# ------------

"""
     residuesdict(residue_list; model=All, chain=All, group=All, residue=All)

This function returns a dictionary (using PDB residue numbers as keys) with the selected
subset of residues. The residues are selected using the keyword arguments `model`, `chain`,
`group` and `residue`. You can use the type `All` (default value) to avoid filtering at
a particular level.
"""
function residuesdict(
    residue_list::AbstractArray{PDBResidue,N};
    model::ResidueQueryTypes = All,
    chain::ResidueQueryTypes = All,
    group::ResidueQueryTypes = All,
    residue::ResidueQueryTypes = All,
) where {N}
    dict = sizehint!(OrderedDict{String,PDBResidue}(), length(residue_list))
    for res in residue_list
        if isresidue(res, model = model, chain = chain, group = group, residue = residue)
            dict[res.id.number] = res
        end
    end
    dict
end

# Deprecation warning for the positional arguments version
function residuesdict(
    residue_list::AbstractArray{PDBResidue,N},
    model,
    chain,
    group,
    residue,
) where {N}
    Base.depwarn(
        "residuesdict using positional arguments is deprecated in favor of keyword arguments: residuesdict(residue_list; model, chain, group, residue)",
        :residuesdict,
        force = true,
    )
    residuesdict(
        residue_list;
        model = model,
        chain = chain,
        group = group,
        residue = residue,
    )
end

"""
`@residuesdict ... model ... chain ... group ... residue ...`

This macro returns a dictionary (using PDB residue numbers as keys) with the selected
subset of residues from a list of residues. You can use the type `All` to avoid filtering
that option.

**DEPRECATED:** This macro is deprecated. Use the [`residuesdict`](@ref) function instead.
"""
macro residuesdict(
    residue_list,
    model::Symbol,
    m,
    chain::Symbol,
    c,
    group::Symbol,
    g,
    residue::Symbol,
    r,
)
    if model == :model && chain == :chain && group == :group && residue == :residue
        Base.depwarn(
            "Using @residuesdict macro is deprecated in favor of residuesdict function with keyword arguments: residuesdict(residue_list; model, chain, group, residue)",
            Symbol("@residuesdict"),
            force = true,
        )
        return :(residuesdict(
            $(esc(residue_list));
            model = $(esc(m)),
            chain = $(esc(c)),
            group = $(esc(g)),
            residue = $(esc(r)),
        ))
    else
        throw(
            ArgumentError(
                "The signature is @residuesdict ___ model ___ chain ___ group ___ residue ___",
            ),
        )
    end
end

# select_atoms
# ------------

"""
    select_atoms(residue_list; model=All, chain=All, group=All, residue=All, atom=All, alt_id=All, charge=All)

This function returns a vector of `PDBAtom`s with the selected subset of atoms from a list
of residues. The atoms are selected using the keyword arguments `model`, `chain`, `group`,
`residue`, `atom`, `alt_id`, and `charge`. You can use the type `All` (default value) to avoid
filtering at a particular level.
"""
function select_atoms(
    residue_list::AbstractArray{PDBResidue,N};
    model::ResidueQueryTypes = All,
    chain::ResidueQueryTypes = All,
    group::ResidueQueryTypes = All,
    residue::ResidueQueryTypes = All,
    atom::ResidueQueryTypes = All,
    alt_id::ResidueQueryTypes = All,
    charge::ResidueQueryTypes = All,
) where {N}
    atom_list = PDBAtom[]
    @inbounds for r in residue_list
        if isresidue(r, model = model, chain = chain, group = group, residue = residue)
            for a in r.atoms
                if isatom(a, atom) && _is(a.alt_id, alt_id) && _is(a.charge, charge)
                    push!(atom_list, a)
                end
            end
        end
    end
    atom_list
end

# Deprecation warning for the positional arguments version
function atoms(
    residue_list::AbstractArray{PDBResidue,N},
    model,
    chain,
    group,
    residue,
    atom,
) where {N}
    Base.depwarn(
        "atoms is deprecated in favor of select_atoms(residue_list; model, chain, group, residue, atom)",
        :atoms,
        force = true,
    )
    select_atoms(
        residue_list;
        model = model,
        chain = chain,
        group = group,
        residue = residue,
        atom = atom,
    )
end

"""
`@atoms ... model ... chain ... group ... residue ... atom ...`

These return a vector of `PDBAtom`s with the selected subset of atoms from a list of
residues. You can use the type `All` to avoid filtering that option.

**DEPRECATED:** This macro is deprecated. Use the [`select_atoms`](@ref) function instead.
"""
macro atoms(
    residue_list,
    model::Symbol,
    m,
    chain::Symbol,
    c,
    group::Symbol,
    g,
    residue::Symbol,
    r,
    atom::Symbol,
    a,
)
    if model == :model &&
       chain == :chain &&
       group == :group &&
       residue == :residue &&
       atom == :atom
        Base.depwarn(
            "Using the @atoms macro is deprecated in favor of the select_atoms function: select_atoms(residue_list; model, chain, group, residue, atom)",
            Symbol("@atoms"),
            force = true,
        )
        return :(select_atoms(
            $(esc(residue_list));
            model = $(esc(m)),
            chain = $(esc(c)),
            group = $(esc(g)),
            residue = $(esc(r)),
            atom = $(esc(a)),
        ))
    else
        throw(
            ArgumentError(
                "The signature is @atoms ___ model ___ chain ___ group ___ residue ___ atom ___",
            ),
        )
    end
end

# Special find...
# ===============

# This _find implementation should be faster than Base.find. It is faster than Base.find
# for small vectors (like atoms) because and allocations is't a problem:
# https://discourse.julialang.org/t/avoid-calling-push-in-base-find/1336
function _find(f::Function, vector::Vector{T}) where {T}
    N = length(vector)
    indices = Array{Int}(undef, N)
    j = 0
    @inbounds for i = 1:N
        if f(vector[i])
            j += 1
            indices[j] = i
        end
    end
    resize!(indices, j)
end

"""
Returns a list with the index of the heavy atoms (all atoms except hydrogen) in
the `PDBResidue`
"""
function findheavy(atoms::Vector{PDBAtom})
    _find(atom -> atom.element != "H", atoms)
end

findheavy(res::PDBResidue) = findheavy(res.atoms)

"""
`findatoms(res::PDBResidue, atom::String)`

Returns a index vector of the atoms with the given `atom` name.
"""
function findatoms(atoms::Vector{PDBAtom}, atom::String)
    _find(a -> a.atom == atom, atoms)
end

findatoms(res::PDBResidue, atom::String) = findatoms(res.atoms, atom)

"""
Returns a vector of indices for `CB` (`CA` for `GLY`)
"""
function findCB(res::PDBResidue)
    atom = res.id.name == "GLY" ? "CA" : "CB"
    findatoms(res, atom)
end

# occupancy
# =========

"""
Takes a `PDBResidue` and a `Vector` of atom indices.
Returns the index value of the `Vector` with maximum occupancy.
"""
function selectbestoccupancy(atoms::Vector{PDBAtom}, indices::Vector{Int})
    Ni = length(indices)
    @assert Ni != 0 "There are no atom indices"
    if Ni == 1
        return (indices[1])
    end
    Na = length(atoms)
    @assert Ni ≤ Na "There are more atom indices ($Ni) than atoms in the Residue ($Na)"
    indice = 0
    occupancy = -Inf
    for i in indices
        actual_occupancy = atoms[i].occupancy
        if actual_occupancy > occupancy
            occupancy = actual_occupancy
            indice = i
        end
    end
    return (indice)
end

selectbestoccupancy(res::PDBResidue, indices) = selectbestoccupancy(res.atoms, indices)

"""
Takes a `Vector` of `PDBAtom`s and returns a `Vector` of the `PDBAtom`s with best occupancy.
"""
function bestoccupancy(atoms::Vector{PDBAtom})::Vector{PDBAtom}
    N = length(atoms)
    if N == 0
        @warn("There are no atoms.")
        return (atoms)
    elseif N == 1
        return (atoms)
    else
        atomdict = sizehint!(OrderedDict{String,PDBAtom}(), N)
        for atom in atoms
            name = atom.atom
            if haskey(atomdict, name)
                if atom.occupancy > atomdict[name].occupancy
                    atomdict[name] = atom
                end
            else
                atomdict[name] = atom
            end
        end
        return (collect(values(atomdict)))
    end
end

function bestoccupancy(res::PDBResidue) # TO DO: Test it!
    new_res = copy(res)
    new_res.atoms = bestoccupancy(res.atoms)
    new_res
end

# More Distances
# ==============

function _update_squared_distance(a_atoms, b_atoms, i, j, dist)
    actual_dist = squared_distance(a_atoms[i], b_atoms[j])
    if actual_dist < dist
        return (actual_dist)
    else
        return (dist)
    end
end

"""
`squared_distance(A::PDBResidue, B::PDBResidue; criteria::String="All")`

Returns the squared distance between the residues `A` and `B`.
The available `criteria` are: `Heavy`, `All`, `CA`, `CB` (`CA` for `GLY`)
"""
function squared_distance(A::PDBResidue, B::PDBResidue; criteria::String = "All")
    a = A.atoms
    b = B.atoms
    dist = Inf
    if criteria == "All"
        Na = length(a)
        Nb = length(b)
        @inbounds for i = 1:Na
            for j = 1:Nb
                dist = _update_squared_distance(a, b, i, j, dist)
            end
        end
    elseif criteria == "Heavy"
        indices_a = findheavy(a)
        indices_b = findheavy(b)
        if length(indices_a) != 0 && length(indices_b) != 0
            for i in indices_a
                for j in indices_b
                    dist = _update_squared_distance(a, b, i, j, dist)
                end
            end
        end
    elseif criteria == "CA"
        indices_a = findatoms(a, "CA")
        indices_b = findatoms(b, "CA")
        if length(indices_a) != 0 && length(indices_b) != 0
            for i in indices_a
                for j in indices_b
                    dist = _update_squared_distance(a, b, i, j, dist)
                end
            end
        end
    elseif criteria == "CB"
        indices_a = findCB(A) # findCB needs residues instead of atoms
        indices_b = findCB(B)
        if length(indices_a) != 0 && length(indices_b) != 0
            for i in indices_a
                for j in indices_b
                    dist = _update_squared_distance(a, b, i, j, dist)
                end
            end
        end
    end
    dist
end

function distance(A::PDBResidue, B::PDBResidue; criteria::String = "All")
    sqrt(squared_distance(A, B, criteria = criteria))
end

"""
`any(f::Function, a::PDBResidue, b::PDBResidue)`

Test if the function `f` is true for any pair of atoms between the residues `a` and `b`
"""
Base.any(f::Function, a::PDBResidue, b::PDBResidue) = any(f, a.atoms, b.atoms)

function Base.any(f::Function, a_atoms::Vector{PDBAtom}, b_atoms::Vector{PDBAtom})
    @inbounds for a in a_atoms
        for b in b_atoms
            if f(a, b)
                return (true)
            end
        end
    end
    return (false)
end

"""
`contact(A::PDBResidue, B::PDBResidue, limit::AbstractFloat; criteria::String="All")`

Returns `true` if the residues `A` and `B` are at contact distance (`limit`).
The available distance `criteria` are: `Heavy`, `All`, `CA`, `CB` (`CA` for `GLY`)
"""
function contact(
    A::PDBResidue,
    B::PDBResidue,
    limit::AbstractFloat;
    criteria::String = "All",
)
    squared_limit = limit^2
    a = A.atoms
    b = B.atoms
    if criteria == "All"
        Na = length(a)
        Nb = length(b)
        @inbounds for i = 1:Na
            ai = a[i]
            for j = 1:Nb
                if _squared_limit_contact(ai, b[j], squared_limit)
                    return (true)
                end
            end
        end
    elseif criteria == "Heavy"
        indices_a = findheavy(a)
        indices_b = findheavy(b)
        if length(indices_a) != 0 && length(indices_b) != 0
            @inbounds for i in indices_a
                ai = a[i]
                for j in indices_b
                    if _squared_limit_contact(ai, b[j], squared_limit)
                        return (true)
                    end
                end
            end
        end
    elseif criteria == "CA"
        indices_a = findatoms(a, "CA")
        indices_b = findatoms(b, "CA")
        if length(indices_a) != 0 && length(indices_b) != 0
            @inbounds for i in indices_a
                ai = a[i]
                for j in indices_b
                    if _squared_limit_contact(ai, b[j], squared_limit)
                        return (true)
                    end
                end
            end
        end
    elseif criteria == "CB"
        indices_a = findCB(A) # findCB needs residues instead of atoms
        indices_b = findCB(B)
        if length(indices_a) != 0 && length(indices_b) != 0
            @inbounds for i in indices_a
                ai = a[i]
                for j in indices_b
                    if _squared_limit_contact(ai, b[j], squared_limit)
                        return (true)
                    end
                end
            end
        end
    end
    false
end

# Vectorize
# =========

# PLM
# ---

"""
It creates a `NamedArray` containing a `PairwiseListMatrix` where each element
(column, row) is identified with a `PDBResidue` from the input vector. You can indicate
the value type of the matrix (default to `Float64`), if the list should have the
diagonal values (default to `Val{false}`) and the diagonal values (default to `NaN`).
"""
function residuepairsmatrix(
    residue_list::Vector{PDBResidue},
    ::Type{T},
    ::Type{Val{diagonal}},
    diagonalvalue::T,
) where {T,diagonal}
    plm = PairwiseListMatrix(T, length(residue_list), diagonal, diagonalvalue)
    resnames = [
        string(res.id.model, '_', res.id.chain, '_', res.id.group, '_', res.id.number)
        for res in residue_list
    ]
    nplm = setlabels(plm, resnames)
    setdimnames!(nplm, ["Res1", "Res2"])
    nplm::NamedArray{
        T,
        2,
        PairwiseListMatrix{T,diagonal,Vector{T}},
        NTuple{2,OrderedDict{String,Int}},
    }
end

function residuepairsmatrix(residue_list::Vector{PDBResidue})
    residuepairsmatrix(residue_list, Float64, Val{false}, NaN)
end

# Contacts and distances
# ----------------------

function squared_distance(residues::Vector{PDBResidue}; criteria::String = "All")
    nplm = residuepairsmatrix(residues, Float64, Val{false}, 0.0)
    plm = getarray(nplm)
    @iterateupper plm false begin
        list[k] = squared_distance(residues[i], residues[j], criteria = criteria)
    end
    nplm
end

"""
`contact(residues::Vector{PDBResidue}, limit::AbstractFloat; criteria::String="All")`

If `contact` takes a `Vector{PDBResidue}`, It returns a matrix with all the pairwise
comparisons (contact map).
"""
function contact(
    residues::Vector{PDBResidue},
    limit::AbstractFloat;
    criteria::String = "All",
)
    nplm = residuepairsmatrix(residues, Bool, Val{false}, true)
    plm = getarray(nplm)
    @iterateupper plm false begin
        list[k] = contact(residues[i], residues[j], limit, criteria = criteria)
    end
    nplm
end

"""
`distance(residues::Vector{PDBResidue}; criteria::String="All")`

If `distance` takes a `Vector{PDBResidue}` returns a `PairwiseListMatrix{Float64, false}`
with all the pairwise comparisons (distance matrix).
"""
function distance(residues::Vector{PDBResidue}; criteria::String = "All")
    nplm = residuepairsmatrix(residues, Float64, Val{false}, 0.0)
    plm = getarray(nplm)
    @iterateupper plm false begin
        list[k] = distance(residues[i], residues[j], criteria = criteria)
    end
    nplm
end

# Proximity average
# -----------------

"""
`proximitymean` calculates the proximity mean/average for each residue as the average
score (from a `scores` list) of all the residues within a certain physical distance to a
given amino acid. The score of that residue is not included in the mean unless you set
`include` to `true`. The default values are 6.05 for the distance threshold/`limit` and
`"Heavy"` for the `criteria` keyword argument. This function allows to calculate pMI
(proximity mutual information) and pC (proximity conservation) as in *Buslje et al.*.

# References

  - [Marino Buslje, Cristina, et al. "Networks of high mutual information define the
    structural proximity of catalytic sites: implications for catalytic residue
    identification." PLoS computational biology 6.11 (2010):
    e1000978.](@cite 10.1371/journal.pcbi.1000978)
"""
function proximitymean(
    residues::Vector{PDBResidue},
    scores::AbstractVector{T},
    limit::T = 6.05;
    criteria::String = "Heavy",
    include::Bool = false,
) where {T<:AbstractFloat}
    N = length(residues)
    @assert N == length(scores) "Vectors must have the same length."
    count = zeros(Int, N)
    sum = zeros(T, N)
    offset = include ? 0 : 1
    @inbounds for i = 1:(N-offset)
        res_i = residues[i]
        for j = (i+offset):N
            if include && (i == j)
                count[i] += 1
                sum[i] += scores[i]
            elseif contact(res_i, residues[j], limit, criteria = criteria)
                count[i] += 1
                count[j] += 1
                sum[i] += scores[j]
                sum[j] += scores[i]
            end
        end
    end
    sum ./ count
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
            if atom.atom in Set{String}(["CE2", "CD2", "CZ2", "CZ3", "CH2", "CE3"])
                push!(plane1, atom)
            end
            if atom.atom in Set{String}(["CG", "CD1", "NE1", "CE2", "CD2"])
                push!(plane2, atom)
            end
        end
        push!(planes, plane1)
        push!(planes, plane2)
    end
    planes
end

function _centre(planes::Vector{Vector{PDBAtom}})
    subset = Vector{PDBAtom}[bestoccupancy(atoms) for atoms in planes]
    polyg =
        Vector{Coordinates}[Coordinates[a.coordinates for a in atoms] for atoms in subset]
    Coordinates[sum(points) ./ Float64(length(points)) for points in polyg]
end

# function _simple_normal_and_centre(atoms::Vector{PDBAtom})
#   atoms = bestoccupancy(atoms)
#   points = Coordinates[ a.coordinates for a in atoms ]
#   (cross(points[2] - points[1], points[3] - points[1]), sum(points)./length(points))
# end

# Show PDB* objects (using Format)
# ====================================

const _Format_ResidueID = FormatExpr("{:>15} {:>15} {:>15} {:>15} {:>15} {:>15}\n")
const _Format_ATOM = FormatExpr("{:>50} {:>15} {:>15} {:>15} {:>15} {:>15} {:>15}\n")

const _Format_ResidueID_Res = FormatExpr("\t\t{:>15} {:>15} {:>15} {:>15} {:>15} {:>15}\n")
const _Format_ATOM_Res =
    FormatExpr("\t\t{:<10} {:>50} {:>15} {:>15} {:>15} {:>15} {:>15} {:>15}\n")

function Base.show(io::IO, id::PDBResidueIdentifier)
    printfmt(
        io,
        _Format_ResidueID,
        "PDBe_number",
        "number",
        "name",
        "group",
        "model",
        "chain",
    )
    printfmt(
        io,
        _Format_ResidueID,
        string('"', id.PDBe_number, '"'),
        string('"', id.number, '"'),
        string('"', id.name, '"'),
        string('"', id.group, '"'),
        string('"', id.model, '"'),
        string('"', id.chain, '"'),
    )
end

function Base.show(io::IO, atom::PDBAtom)
    printfmt(
        io,
        _Format_ATOM,
        "coordinates",
        "atom",
        "element",
        "occupancy",
        "B",
        "alt_id",
        "charge",
    )
    printfmt(
        io,
        _Format_ATOM,
        atom.coordinates,
        string('"', atom.atom, '"'),
        string('"', atom.element, '"'),
        atom.occupancy,
        string('"', atom.B, '"'),
        string('"', atom.alt_id, '"'),
        string('"', atom.charge, '"'),
    )
end

function Base.show(io::IO, res::PDBResidue)
    println(io, "PDBResidue:\n\tid::PDBResidueIdentifier")
    printfmt(
        io,
        _Format_ResidueID_Res,
        "PDBe_number",
        "number",
        "name",
        "group",
        "model",
        "chain",
    )
    printfmt(
        io,
        _Format_ResidueID_Res,
        string('"', res.id.PDBe_number, '"'),
        string('"', res.id.number, '"'),
        string('"', res.id.name, '"'),
        string('"', res.id.group, '"'),
        string('"', res.id.model, '"'),
        string('"', res.id.chain, '"'),
    )
    len = length(res)
    println(io, "\tatoms::Vector{PDBAtom}\tlength: ", len)
    for i = 1:len
        printfmt(
            io,
            _Format_ATOM_Res,
            "",
            "coordinates",
            "atom",
            "element",
            "occupancy",
            "B",
            "alt_id",
            "charge",
        )
        printfmt(
            io,
            _Format_ATOM_Res,
            string(i, ":"),
            res.atoms[i].coordinates,
            string('"', res.atoms[i].atom, '"'),
            string('"', res.atoms[i].element, '"'),
            res.atoms[i].occupancy,
            string('"', res.atoms[i].B, '"'),
            string('"', res.atoms[i].alt_id, '"'),
            string('"', res.atoms[i].charge, '"'),
        )
    end
end
