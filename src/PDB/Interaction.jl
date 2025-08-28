_with_vdw(a::PDBAtom, resname_a::String) = (resname_a, a.atom) in keys(vanderwaalsradius)

_with_cov(a::PDBAtom, resname_a::String) = a.element in keys(covalentradius)

# Detect peptide bonds between two atoms
"""
    peptide_bond(
        res_a::PDBResidue,
        a::PDBAtom,
        res_b::PDBResidue,
        b::PDBAtom,
        cutoff = 1.38,
    )

Return `true` if the atoms form a peptide bond.

The bond is recognized when the atoms are named `"C"` and `"N"`, belong to the
same chain and model, and are separated by at most `cutoff` Å. If neither atom is
`"C"` nor `"N"`, the function returns `missing`.

Peptide bonds have been reported in the 1.28–1.38 Å range, so the maximum value
is used to be as inclusive as possible.

# References

  - [Panjikar, Santosh, and Manfred S. Weiss. "Peptide bonds revisited." IUCrJ
    12.3 (2025): 307–321.](@cite 10.1107/S2052252525002106)
"""
function peptide_bond(
    res_a::PDBResidue,
    a::PDBAtom,
    res_b::PDBResidue,
    b::PDBAtom,
    cutoff = 1.38,
)
    names = (a.atom, b.atom)
    if !("C" in names && "N" in names)
        return missing
    end
    if res_a.id.chain == res_b.id.chain &&
       res_a.id.model == res_b.id.model &&
       distance(a, b) <= cutoff
        return true
    end
    return false
end

peptide_bond(res_a::PDBResidue, ia::Int, res_b::PDBResidue, ib::Int, cutoff = 1.38) =
    peptide_bond(res_a, res_a.atoms[ia], res_b, res_b.atoms[ib], cutoff)

ishydrophobic(a::PDBAtom, resname_a::String) = (resname_a, a.atom) in _hydrophobic

"""
Returns true if the atom, e.g. `("HIS","CG")`, is an aromatic atom in the residue.
"""
isaromatic(a::PDBAtom, resname_a::String) = (resname_a, a.atom) in _aromatic

"""
Returns true if the atom, e.g. `("ARG","NE")`, is a cationic atom in the residue.
"""
iscationic(a::PDBAtom, resname_a::String) = (resname_a, a.atom) in _cationic

"""
Returns true if the atom, e.g. `("GLU","CD")`, is an anionic atom in the residue.
"""
isanionic(a::PDBAtom, resname_a::String) = (resname_a, a.atom) in _anionic

"""
Returns true if the atom, e.g. `("ARG","N")`, is a donor in H bonds.
"""
ishbonddonor(a::PDBAtom, resname_a::String) = (resname_a, a.atom) in keys(_hbond_donor)

"""
Returns true if the atom, e.g. `("ARG","O")`, is an acceptor in H bonds.
"""
ishbondacceptor(a::PDBAtom, resname_a::String) =
    (resname_a, a.atom) in keys(_hbond_acceptor)

"""
`any(f::Function, a::PDBResidue, b::PDBResidue, criteria::Function)`

Test if the function `f` is true for any pair of atoms between the residues `a` and `b`.
This function only test atoms that returns `true` for the fuction `criteria`.
"""
function Base.any(f::Function, a::PDBResidue, b::PDBResidue, criteria::Function)
    resname_a, resname_b = a.id.name, b.id.name
    a_atoms = a.atoms
    b_atoms = b.atoms
    indices_a = _find(x -> criteria(x, resname_a), a_atoms)
    indices_b = _find(x -> criteria(x, resname_b), b_atoms)
    if length(indices_a) != 0 && length(indices_b) != 0
        @inbounds for i in indices_a
            for j in indices_b
                if f(a_atoms[i], b_atoms[j], resname_a, resname_b)
                    return (true)
                end
            end
        end
    end
    return (false)
end

# Interaction types
# =================

# van der Waals
# -------------

"""
Test if two atoms or residues are in van der Waals contact using:
`distance(a,b) <= 0.5 + vanderwaalsradius[a] + vanderwaalsradius[b]`.
It returns distance `<= 0.5` if the atoms aren't in `vanderwaalsradius`.
"""
function vanderwaals(a::PDBAtom, b::PDBAtom, resname_a, resname_b)
    return (
        distance(a, b) <=
        0.5 +
        get(vanderwaalsradius, (resname_a, a.atom), 0.0) +
        get(vanderwaalsradius, (resname_b, b.atom), 0.0)
    )
end

vanderwaals(a::PDBResidue, b::PDBResidue) = any(vanderwaals, a, b, _with_vdw)

# van der Waals clash
# -------------------

"""
    vanderwaalsclash(res_a::PDBResidue, a::PDBAtom, res_b::PDBResidue, b::PDBAtom)

Return `true` if the distance between the atoms is less than the sum of their
`vanderwaalsradius` values.

Pairs detected as peptide bonds by [`peptide_bond`](@ref) are ignored. If the
atoms are not listed (for example `OXT`), the radius of the element is used.
Unknown elements fall back to `0.0`. Only distances are checked; no chemical
context is considered.
"""

function vanderwaalsclash(res_a::PDBResidue, a::PDBAtom, res_b::PDBResidue, b::PDBAtom)
    bond = peptide_bond(res_a, a, res_b, b)
    bond === true && return false
    resname_a = res_a.id.name
    resname_b = res_b.id.name
    return (
        distance(a, b) <=
        get(
            vanderwaalsradius,
            (resname_a, a.atom),
            get(vanderwaalsradius, (resname_a, a.element), 0.0),
        ) + get(
            vanderwaalsradius,
            (resname_b, b.atom),
            get(vanderwaalsradius, (resname_b, b.element), 0.0),
        )
    )
end

vanderwaalsclash(res_a::PDBResidue, ia::Int, res_b::PDBResidue, ib::Int) =
    vanderwaalsclash(res_a, res_a.atoms[ia], res_b, res_b.atoms[ib])

function vanderwaalsclash(a::PDBAtom, b::PDBAtom, resname_a, resname_b)
    Base.depwarn(
        "vanderwaalsclash(a,b,resname_a,resname_b) is deprecated. " *
        "Use vanderwaalsclash(a,b,residue_a,residue_b) instead.",
        :vanderwaalsclash,
    )
    return (
        distance(a, b) <=
        get(
            vanderwaalsradius,
            (resname_a, a.atom),
            get(vanderwaalsradius, (resname_a, a.element), 0.0),
        ) + get(
            vanderwaalsradius,
            (resname_b, b.atom),
            get(vanderwaalsradius, (resname_b, b.element), 0.0),
        )
    )
end

function vanderwaalsclash(a::PDBResidue, b::PDBResidue)
    resname_a, resname_b = a.id.name, b.id.name
    a_atoms = a.atoms
    b_atoms = b.atoms
    indices_a = _find(x -> _with_vdw(x, resname_a), a_atoms)
    indices_b = _find(x -> _with_vdw(x, resname_b), b_atoms)
    if !isempty(indices_a) && !isempty(indices_b)
        @inbounds for i in indices_a
            for j in indices_b
                atom_a = a_atoms[i]
                atom_b = b_atoms[j]
                bond = peptide_bond(a, atom_a, b, atom_b)
                if (bond !== true) && vanderwaalsclash(a, atom_a, b, atom_b)
                    return true
                end
            end
        end
    end
    return false
end

# Covalent
# --------

_covalent(a::PDBAtom, b::PDBAtom, resname_a::String, resname_b::String) = (
    distance(a, b) <=
    get(covalentradius, a.element, 0.0) + get(covalentradius, b.element, 0.0)
)

function covalent(a::PDBAtom, b::PDBAtom, resname_a, resname_b)
    Base.depwarn(
        "covalent(a,b,resname_a,resname_b) is deprecated. " *
        "Use covalent(res_a,a,res_b,b) instead.",
        :covalent,
    )
    return _covalent(a, b, resname_a, resname_b)
end


covalent(a::PDBAtom, b::PDBAtom) = _covalent(a, b, "", "")

"""
    covalent(res_a::PDBResidue, a::PDBAtom, res_b::PDBResidue, b::PDBAtom)

Return `true` if the distance between atoms is less than the sum of the
`covalentradius` of each atom.

For residues, the check iterates over atoms with known covalent radii and also
reports `true` when a peptide bond is detected by [`peptide_bond`](@ref).
Peptide bonds are considered covalent even when the `C` and `N` atoms are
farther apart than the sum of their radii.

!!! warning

    Atom pairs separated by less than the sum of their covalent radii are
    reported as covalent, even if they correspond to steric clashes rather than
    true bonds.

This method only verifies distances and does not inspect bonding angles or other
chemical context.
"""
function covalent(res_a::PDBResidue, a::PDBAtom, res_b::PDBResidue, b::PDBAtom)
    bond = peptide_bond(res_a, a, res_b, b)
    bond === true && return true
    return covalent(a, b, res_a.id.name, res_b.id.name)
end

covalent(res_a::PDBResidue, ia::Int, res_b::PDBResidue, ib::Int) =
    covalent(res_a, res_a.atoms[ia], res_b, res_b.atoms[ib])

function covalent(a::PDBResidue, b::PDBResidue)
    a_atoms = a.atoms
    b_atoms = b.atoms
    indices_a = _find(x -> _with_cov(x, a.id.name), a_atoms)
    indices_b = _find(x -> _with_cov(x, b.id.name), b_atoms)
    if !isempty(indices_a) && !isempty(indices_b)
        @inbounds for i in indices_a
            for j in indices_b
                if covalent(a, a_atoms[i], b, b_atoms[j])
                    return true
                end
            end
        end
    end
    return false
end

# Disulphide
# ----------

_issulphurcys(a::PDBAtom, resname_a) = resname_a == "CYS" && a.element == "S"

"""
Returns `true` if two `CYS`'s `S` are at 2.08 Å or less
"""
function disulphide(a::PDBAtom, b::PDBAtom, resname_a, resname_b)
    if _issulphurcys(a, resname_a) && _issulphurcys(b, resname_b)
        return (squared_distance(a, b) <= (2.08^2))
    end
    return (false)
end

disulphide(a::PDBResidue, b::PDBResidue) = any(disulphide, a, b, _issulphurcys)

# Aromatic-Sulphur
# ----------------

_issulphur(a::PDBAtom) = a.element == "S"

"""
Returns `true` if an sulphur and an aromatic atoms are 5.3 Å or less"
"""
function aromaticsulphur(a::PDBAtom, b::PDBAtom, resname_a, resname_b)
    if (_issulphur(a) && isaromatic(b, resname_b)) ||
       (_issulphur(b) && isaromatic(a, resname_a))
        return (squared_distance(a, b) <= 28.09) # 28.09 == 5.3 ^ 2
    end
    return (false)
end

_issulphuroraromatic(a::PDBAtom, resname_a) = _issulphur(a) || isaromatic(a, resname_a)

aromaticsulphur(a::PDBResidue, b::PDBResidue) =
    any(aromaticsulphur, a, b, _issulphuroraromatic)

# Π-Cation
# --------

"""
There's a Π-Cation interaction if a cationic and an aromatic atoms are at 6.0 Å or less
"""
function pication(a::PDBAtom, b::PDBAtom, resname_a, resname_b)
    if (iscationic(a, resname_a) && isaromatic(b, resname_b)) ||
       (iscationic(b, resname_b) && isaromatic(a, resname_a))
        return (squared_distance(a, b) <= 36.0) # 36.0 == 6.0 ^ 2
    end
    return (false)
end

_iscationicoraromatic(a::PDBAtom, resname_a) =
    iscationic(a, resname_a) || isaromatic(a, resname_a)

pication(a::PDBResidue, b::PDBResidue) = any(pication, a, b, _iscationicoraromatic)

# Aromatic
# --------

"""
There's an aromatic interaction if centriods are at 6.0 Å or less.
"""
function aromatic(a::PDBResidue, b::PDBResidue)
    threshold = 36.0 # 6.0 ^ 2
    if (a.id.name in _aromatic_res) &&
       (b.id.name in _aromatic_res) &&
       (squared_distance(a, b) <= threshold)
        centres_a = _centre(_get_plane(a))
        centres_b = _centre(_get_plane(b))
        return (any(
            squared_distance(centroid_a, centroid_b) <= threshold for
            centroid_a in centres_a, centroid_b in centres_b
        ))
    end
    return (false)
end

# Ionic
# -----

"""
There's an ionic interaction if a cationic and an anionic atoms are at 6.0 Å or less.
"""
function ionic(a::PDBAtom, b::PDBAtom, resname_a, resname_b)
    if (iscationic(a, resname_a) && isanionic(b, resname_b)) ||
       (iscationic(b, resname_b) && isanionic(a, resname_a))
        return (squared_distance(a, b) <= 36.0) # 36.0 == 6.0 ^ 2
    end
    return (false)
end

_iscationicoranionic(a::PDBAtom, resname_a) =
    iscationic(a, resname_a) || isanionic(a, resname_a)

ionic(a::PDBResidue, b::PDBResidue) = any(ionic, a, b, _iscationicoranionic)

# Hydrophobic contact
# -------------------

"""
There's an hydrophobic interaction if two hydrophobic atoms are at 5.0 Å or less.
"""
function hydrophobic(a::PDBAtom, b::PDBAtom, resname_a, resname_b)
    if ishydrophobic(a, resname_a) && ishydrophobic(b, resname_b)
        return (squared_distance(a, b) <= 25.0) # 5.0 ^ 2
    end
    return (false)
end

hydrophobic(a::PDBResidue, b::PDBResidue) = any(hydrophobic, a, b, ishydrophobic)

# Hydrogen bonds
# --------------

function _find_antecedent(res::PDBResidue, a::PDBAtom)
    ids = _hbond_acceptor[(res.id.name, a.atom)]
    _find(a -> a.atom in ids, res.atoms)
end

function _find_h(res::PDBResidue, a::PDBAtom)
    ids = _hbond_donor[(res.id.name, a.atom)]
    _find(a -> a.atom in ids, res.atoms)
end

function _hbond_kernel(donor, acceptor, indices_donor, indices_acceptor)
    @inbounds for i in indices_donor
        don = donor.atoms[i]
        indices_h = _find_h(donor, don)
        if length(indices_h) == 0
            continue
        end
        for j in indices_acceptor
            acc = acceptor.atoms[j]
            indices_ant = _find_antecedent(acceptor, acc)
            if squared_distance(don, acc) <= (3.9^2) && length(indices_ant) != 0
                for k in indices_h
                    hyd = donor.atoms[k]
                    if squared_distance(hyd, acc) <= 6.25 && angle(don, hyd, acc) >= 90.0 # 6.25 == 2.5²
                        for ant in indices_ant
                            if angle(don, acc, acceptor.atoms[ant]) >= 90.0 &&
                               angle(hyd, acc, acceptor.atoms[ant]) >= 90.0
                                return (true)
                            end
                        end
                    end
                end
            end
        end
    end
    return (false)
end

function _hydrogenbond_don_acc(donor::PDBResidue, acceptor::PDBResidue)
    if donor != acceptor
        indices_donor = findall(x -> ishbonddonor(x, donor.id.name), donor.atoms)
        indices_acceptor =
            findall(x -> ishbondacceptor(x, acceptor.id.name), acceptor.atoms)
        if length(indices_donor) != 0 && length(indices_acceptor) != 0
            return (_hbond_kernel(donor, acceptor, indices_donor, indices_acceptor))
        end
    end
    return (false)
end


"""
This function only works if there are hydrogens in the structure.
The criteria for a hydrogen bond are:

  - d(Ai, Aj) < 3.9Å
  - d(Ah, Aacc) < 2.5Å
  - θ(Adon, Ah, Aacc) > 90°
  - θ(Adon, Aacc, Aacc-antecedent) > 90°
  - θ(Ah, Aacc, Aacc-antecedent) > 90°

Where Ah is the donated hydrogen atom, Adon is the hydrogen bond donor atom,
Aacc is the hydrogen bond acceptor atom and Aacc-antecednt is the atom antecedent to the
hydrogen bond acceptor atom.
"""
hydrogenbond(a::PDBResidue, b::PDBResidue) =
    _hydrogenbond_don_acc(a, b) || _hydrogenbond_don_acc(b, a)
