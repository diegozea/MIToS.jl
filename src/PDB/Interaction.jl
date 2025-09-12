_with_vdw(a::PDBAtom, resname_a::String) = a.element in keys(VAN_DER_WAALS_RADII)

_with_cov(a::PDBAtom, resname_a::String) = a.element in keys(COVALENT_RADII)

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
    # NOTE: 1.38 is lower than the bond length calculated using covalent radii, 
    # which is 1.617 Å calculated as (0.76 Å + 0.71 Å) * 1.1
    # This means, that the covalent function will always return true for peptide bonds.
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

function peptide_bond(res_a::PDBResidue, res_b::PDBResidue, cutoff = 1.38)
    if res_a != res_b
        indices_a = findall(x -> x.atom == "C", res_a.atoms)
        indices_b = findall(x -> x.atom == "N", res_b.atoms)
        if length(indices_a) != 0 && length(indices_b) != 0
            @inbounds for i in indices_a
                for j in indices_b
                    if peptide_bond(res_a, res_a.atoms[i], res_b, res_b.atoms[j], cutoff)
                        return true
                    end
                end
            end
            return false
        end
    end
    missing
end

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
`any(f::Function, a::PDBResidue, b::PDBResidue, criteria::Function; kwargs...)`

Test if the function `f` is true for any pair of atoms between the residues `a` and `b`.
This function only tests atoms that return `true` for the function `criteria`.
Any keyword arguments provided are forwarded to `f`. The function `f` must return a
boolean value and have the following signature:

```julia
f(a::PDBAtom, b::PDBAtom, resname_a::String, resname_b::String; kwargs...)
```

For example, the functions [`covalent`](@ref), [`vanderwaals`](@ref),
[`vanderwaalsclash`](@ref), [`disulphide`](@ref), [`aromaticsulphur`](@ref),
[`pication`](@ref), [`ionic`](@ref), and [`hydrophobic`](@ref) can be used with `any`.

The function `criteria` should also return a boolean value and it should have the
following signature:

```julia
criteria(a::PDBAtom, resname_a::String)
```
"""
function Base.any(f::Function, a::PDBResidue, b::PDBResidue, criteria::Function; kwargs...)
    resname_a, resname_b = a.id.name, b.id.name
    a_atoms = a.atoms
    b_atoms = b.atoms
    indices_a = _find(x -> criteria(x, resname_a), a_atoms)
    indices_b = _find(x -> criteria(x, resname_b), b_atoms)
    if length(indices_a) != 0 && length(indices_b) != 0
        @inbounds for i in indices_a
            for j in indices_b
                if f(a_atoms[i], b_atoms[j], resname_a, resname_b; kwargs...)
                    return true
                end
            end
        end
    end
    false
end

# Interaction types
# =================

# Covalent
# --------

function _get_covalent_radius(atom::PDBAtom)
    element = atom.element
    if haskey(COVALENT_RADII, element)
        COVALENT_RADII[element]
    else
        @warn(
            "Element $element not found in `COVALENT_RADII`; using 0.0 as default",
            maxlog = 1,
            _id = ("COVALENT_RADII", element)
        )
        0.0
    end
end

"""
    covalent(
        a::PDBAtom,
        b::PDBAtom;
        tolerance_factor::Float64 = 1.1,
    )

Return `true` if the distance between atoms is less than the sum of the
[`COVALENT_RADII`](@ref) of each atom taking into account the tolerance factor.
By default, the tolerance factor is set to `1.1`, allowing for a 10% increase in the sum
of the covalent radii. This multiplicative factor can be adjusted using the
`tolerance_factor` keyword argument. This function is based on equation 5 from
*Kim and Kim (2015)*. Covalent radii are obtained from Cordero et al. (2008). If the
element is not listed in [`COVALENT_RADII`](@ref), the radius is set to `0.0`, and this
function will return `false`.

!!! warning

    This check considers only interatomic distances; it does not evaluate bond angles,
    coordination, or other chemical context. Any atom pair closer than the sum of their
    covalent radii (taking into account the tolerance factor) is flagged as a covalent
    contact, even if this corresponds to a steric clash rather than a true bond.

# References

    - [Kim, Y. and Kim, W.Y. (2015), Universal Structure Conversion Method for Organic
      Molecules: From Atomic Connectivity to Three-Dimensional Geometry. 
      Bull. Korean Chem. Soc., 36: 1769-1777.](@cite https://doi.org/10.1002/bkcs.10334)
"""
function covalent(
    a::PDBAtom,
    b::PDBAtom;
    tolerance_factor::Float64 = 1.1, # 10% tolerance
)
    # 0.0 as default to avoid defining bonds for elements without known covalent radius    
    rᵢ = _get_covalent_radius(a)
    rⱼ = _get_covalent_radius(b)
    if rᵢ == 0.0 || rⱼ == 0.0
        @warn(
            "Covalent bonds are undefined when the covalent radius is unknown (i.e., 0.0).",
            maxlog = 1
        )
        return false
    end
    dᵢⱼ = distance(a, b)
    dᵢⱼ <= (tolerance_factor * (rᵢ + rⱼ))
end

function covalent(
    a::PDBAtom,
    b::PDBAtom,
    # the `any` function will call this function with the resnames
    resname_a::String,
    resname_b::String;
    tolerance_factor::Float64 = 1.1,
)
    covalent(a, b, tolerance_factor = tolerance_factor)
end

covalent(
    res_a::PDBResidue,
    ia::Int,
    res_b::PDBResidue,
    ib::Int;
    tolerance_factor::Float64 = 1.1,
) = covalent(res_a.atoms[ia], res_b.atoms[ib], tolerance_factor = tolerance_factor)

covalent(
    res_a::PDBResidue,
    atom_a::PDBAtom,
    res_b::PDBResidue,
    atom_b::PDBAtom;
    tolerance_factor::Float64 = 1.1,
) = covalent(atom_a, atom_b, tolerance_factor = tolerance_factor)

function covalent(a::PDBResidue, b::PDBResidue; tolerance_factor::Float64 = 1.1)
    any(covalent, a, b, _with_cov; tolerance_factor = tolerance_factor)
end

# van der Waals
# -------------

function _get_vanderwaals_radius(atom::PDBAtom)
    if haskey(VAN_DER_WAALS_RADII, atom.element)
        VAN_DER_WAALS_RADII[atom.element]
    else
        @warn(
            "Element $(atom.element) not found in `VAN_DER_WAALS_RADII`; using 0.0 as default",
            maxlog = 1,
            _id = ("VAN_DER_WAALS_RADII", atom.element)
        )
        0.0
    end
end

"""
    vanderwaals(a::PDBAtom, b::PDBAtom)

Test if two atoms or residues are in van der Waals contact using the criterion in
Alvarez (2013). That means, if the distance between two atoms is between 0.7 Å less
and 0.7 Å more than the sum of their van der Waals radii. If the element is not listed
in [`VAN_DER_WAALS_RADII`](@ref), the radius is set to `0.0`, and this function will
return `false`.

# References

    - [Alvarez, Santiago. "A cartography of the van der Waals territories." Dalton 
      Transactions 42.24 (2013): 8617-8636.](@cite C3DT50599E)
"""
function vanderwaals(a::PDBAtom, b::PDBAtom)
    d_ab = distance(a, b)
    r_vdW_a = _get_vanderwaals_radius(a)
    r_vdW_b = _get_vanderwaals_radius(b)
    if r_vdW_a == 0.0 || r_vdW_b == 0.0
        false
    else
        dw = r_vdW_a + r_vdW_b
        (dw - 0.7) <= d_ab <= (dw + 0.7)
    end
end

# Definition to use with the `any` function
vanderwaals(a::PDBAtom, b::PDBAtom, resname_a::String, resname_b::String) =
    vanderwaals(a, b)

vanderwaals(a::PDBResidue, b::PDBResidue) = any(vanderwaals, a, b, _with_vdw)

# van der Waals clash
# -------------------

"""
    vanderwaalsclash(a::PDBAtom, b::PDBAtom; tolerance_value::Float64 = -0.7)

This function considered a van der Waals clash if the distance between two atoms falls
within the van der Waals gap as defined by Alvarez (2013). That means, if the distance
between two atoms that are not covalently bonded is less than or equal to the sum of their
van der Waals radii less 0.7 Å (as the default value for `tolerance_value` is `-0.7`).
Here we use the [`covalent`](@ref) function to check for covalent bonds instead of fixing
a distance threshold as in Alvarez (2013) — that paper suggest that _"distances shorter
than the radii sum by more than 1.3 Å correspond most likely to a chemical bond"_. Unknown
elements fall back to `0.0` returning `false`. Only distances are checked; no chemical
context is considered.

!!! warning

    The `tolerance_value` as been set to `-0.7` by default in MIToS 3.2.0 and later.
    Previous versions used `0.0` as default, matching the definition on
    Bickerton et al. (2012). However, that value is too high and results in many
    van der Waals contacts being classified as clashes. MIToS 3.2.0 also replaced the
    van der Waals radii from *Bickerton et al.* (2012) with the values reported by
    *Alvarez* (2013). Alvarez (2013) suggests that _"the position of the shortest
    non-bonded distance can be roughly estimated to be 0.7 Å shorter than the radii
    sum."_

# References

    - [Alvarez, Santiago. "A cartography of the van der Waals territories." Dalton 
      Transactions 42.24 (2013): 8617-8636.](@cite C3DT50599E)
    - [Bickerton, George R., Alicia P. Higueruelo, and Tom L. Blundell. "Comprehensive, 
      atomic-level characterization of structurally characterized protein-protein 
      interactions: the PICCOLO database." BMC bioinformatics 
      12 (2011): 1-15.](@cite 10.1186/1471-2105-12-313)
"""
function vanderwaalsclash(a::PDBAtom, b::PDBAtom; tolerance_value::Float64 = -0.7)
    if covalent(a, b)
        return false
    end
    r_vdW_a = _get_vanderwaals_radius(a)
    r_vdW_b = _get_vanderwaals_radius(b)
    if r_vdW_a == 0.0 || r_vdW_b == 0.0
        false
    else
        distance(a, b) < (r_vdW_a + r_vdW_b + tolerance_value)
    end
end

vanderwaalsclash(
    a::PDBAtom,
    b::PDBAtom,
    resname_a::String,
    resname_b::String;
    tolerance_value::Float64 = -0.7,
) = vanderwaalsclash(a, b; tolerance_value = tolerance_value)

vanderwaalsclash(
    res_a::PDBResidue,
    ia::Int,
    res_b::PDBResidue,
    ib::Int;
    tolerance_value::Float64 = -0.7,
) = vanderwaalsclash(res_a.atoms[ia], res_b.atoms[ib]; tolerance_value = tolerance_value)

vanderwaalsclash(
    res_a::PDBResidue,
    atom_a::PDBAtom,
    res_b::PDBResidue,
    atom_b::PDBAtom;
    tolerance_value::Float64 = -0.7,
) = vanderwaalsclash(atom_a, atom_b; tolerance_value = tolerance_value)

function vanderwaalsclash(a::PDBResidue, b::PDBResidue; tolerance_value::Float64 = -0.7)
    any(vanderwaalsclash, a, b, _with_vdw; tolerance_value = tolerance_value)
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
