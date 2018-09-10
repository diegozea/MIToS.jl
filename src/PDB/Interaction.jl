_with_vdw(a::PDBAtom, resname_a::String) = (resname_a, a.atom) in keys(vanderwaalsradius)

_with_cov(a::PDBAtom, resname_a::String) = a.element in keys(covalentradius)

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
ishbondacceptor(a::PDBAtom, resname_a::String) = (resname_a, a.atom) in keys(_hbond_acceptor)

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
                    return(true)
                end
            end
        end
    end
    return(false)
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
    return( distance(a,b) <= 0.5 +
               get(vanderwaalsradius, (resname_a, a.atom), 0.0) +
               get(vanderwaalsradius, (resname_b, b.atom), 0.0) )
end

vanderwaals(a::PDBResidue, b::PDBResidue) = any(vanderwaals, a, b, _with_vdw)

# van der Waals clash
# -------------------

"""
Returns `true` if the distance between the atoms is less than the sum of the
`vanderwaalsradius` of the atoms. If the atoms aren't on the list (i.e. `OXT`), the
`vanderwaalsradius` of the element is used. If there is not data in the dict,
distance `0.0` is used.
"""
function vanderwaalsclash(a::PDBAtom, b::PDBAtom, resname_a, resname_b)
    return( distance(a,b) <= get(vanderwaalsradius, (resname_a, a.atom),
                                 get(vanderwaalsradius, (resname_a, a.element), 0.0)) +
                             get(vanderwaalsradius, (resname_b, b.atom),
                                 get(vanderwaalsradius, (resname_b, b.element), 0.0)) )
end

vanderwaalsclash(a::PDBResidue, b::PDBResidue) = any(vanderwaalsclash, a, b, _with_vdw)

# Covalent
# --------

"""
Returns `true` if the distance between atoms is less than the sum of the `covalentradius`
of each atom.
"""
function covalent(a::PDBAtom, b::PDBAtom, resname_a, resname_b) # any(... calls it with the res names
    return( distance(a,b) <= get(covalentradius, a.element, 0.0) +
               get(covalentradius, b.element, 0.0) )
end

covalent(a::PDBAtom, b::PDBAtom) = covalent(a, b, "", "")

covalent(a::PDBResidue, b::PDBResidue) = any(covalent, a, b, _with_cov)

# Disulphide
# ----------

_issulphurcys(a::PDBAtom, resname_a) = resname_a == "CYS" && a.element == "S"

"Returns `true` if two `CYS`'s `S` are at 2.08 Å or less"
function disulphide(a::PDBAtom, b::PDBAtom, resname_a, resname_b)
    if _issulphurcys(a, resname_a) && _issulphurcys(b, resname_b)
        return(squared_distance(a,b) <= (2.08 ^ 2))
    end
    return(false)
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
        return(squared_distance(a,b) <= 28.09) # 28.09 == 5.3 ^ 2
    end
    return(false)
end

_issulphuroraromatic(a::PDBAtom, resname_a) = _issulphur(a) || isaromatic(a, resname_a)

aromaticsulphur(a::PDBResidue, b::PDBResidue) = any(aromaticsulphur, a, b, _issulphuroraromatic)

# Π-Cation
# --------

"""
There's a Π-Cation interaction if a cationic and an aromatic atoms are at 6.0 Å or less
"""
function pication(a::PDBAtom, b::PDBAtom, resname_a, resname_b)
    if (iscationic(a, resname_a) && isaromatic(b, resname_b)) ||
            (iscationic(b, resname_b) && isaromatic(a, resname_a))
        return(squared_distance(a,b) <= 36.0) # 36.0 == 6.0 ^ 2
    end
    return(false)
end

_iscationicoraromatic(a::PDBAtom, resname_a) = iscationic(a, resname_a) || isaromatic(a, resname_a)

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
        return(any(squared_distance(centroid_a,centroid_b) <= threshold for centroid_a in centres_a, centroid_b in centres_b))
    end
    return(false)
end

# Ionic
# -----

"""
There's an ionic interaction if a cationic and an anionic atoms are at 6.0 Å or less.
"""
function ionic(a::PDBAtom, b::PDBAtom, resname_a, resname_b)
    if (iscationic(a, resname_a) && isanionic(b, resname_b)) ||
            (iscationic(b, resname_b) && isanionic(a, resname_a))
        return(squared_distance(a,b) <= 36.0) # 36.0 == 6.0 ^ 2
    end
    return(false)
end

_iscationicoranionic(a::PDBAtom, resname_a) = iscationic(a, resname_a) || isanionic(a, resname_a)

ionic(a::PDBResidue, b::PDBResidue) = any(ionic, a, b, _iscationicoranionic)

# Hydrophobic contact
# -------------------

"""
There's an hydrophobic interaction if two hydrophobic atoms are at 5.0 Å or less.
"""
function hydrophobic(a::PDBAtom, b::PDBAtom, resname_a, resname_b)
    if ishydrophobic(a, resname_a) && ishydrophobic(b, resname_b)
        return(squared_distance(a,b) <= 25.0) # 5.0 ^ 2
    end
    return(false)
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
                            if angle(don, acc, acceptor.atoms[ant]) >= 90.0 && angle(hyd, acc, acceptor.atoms[ant]) >= 90.0
                                return(true)
                            end
                        end
                    end
                end
            end
        end
    end
    return(false)
end

function _hydrogenbond_don_acc(donor::PDBResidue, acceptor::PDBResidue)
    if donor != acceptor
        indices_donor = findall(x -> ishbonddonor(x, donor.id.name), donor.atoms)
        indices_acceptor = findall(x -> ishbondacceptor(x, acceptor.id.name), acceptor.atoms)
        if length(indices_donor) != 0 && length(indices_acceptor) != 0
            return(_hbond_kernel(donor, acceptor, indices_donor, indices_acceptor))
        end
    end
    return(false)
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
hydrogenbond(a::PDBResidue, b::PDBResidue) = _hydrogenbond_don_acc(a,b) || _hydrogenbond_don_acc(b,a)
