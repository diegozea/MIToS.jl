import Base: any

__with_vdw(a::PDBAtom) = (a.residueid.name, a.atomid) in keys(vanderwaalsradius)

__with_cov(a::PDBAtom) = a.element in keys(covalentradius)

ishydrophobic(a::PDBAtom) = (a.residueid.name, a.atomid) in __hydrophobic

isaromatic(a::PDBAtom) = (a.residueid.name, a.atomid) in __aromatic

iscationic(a::PDBAtom) = (a.residueid.name, a.atomid) in __cationic

isanionic(a::PDBAtom) = (a.residueid.name, a.atomid) in __anionic

ishbonddonor(a::PDBAtom) = (a.residueid.name, a.atomid) in __hbond_donor

ishbondacceptor(a::PDBAtom) = (a.residueid.name, a.atomid) in __hbond_acceptor

"""Test if the function f is true for any pair of atoms between the residues a and b,
only test atoms that returns true for the fuction criteria"""
function any(f::Function, a::PDBResidue, b::PDBResidue, criteria::Function)
  indices_a = find(criteria, a.atoms)
  indices_b = find(criteria, b.atoms)
  if length(indices_a) != 0 && length(indices_b) != 0
    @inbounds for i in indices_a
      for j in indices_b
        if f(a.atoms[i], b.atoms[j])
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

"""Returns dist <= 0.5 if the atoms aren't in vanderwaalsradius"""
function vanderwaals(a::PDBAtom, b::PDBAtom)
  return( distance(a,b) <= 0.5 +
          get(vanderwaalsradius, (a.residueid.name, a.atomid), 0.0) +
          get(vanderwaalsradius, (b.residueid.name, b.atomid), 0.0) )
end

vanderwaals(a::PDBResidue, b::PDBResidue) = any(vanderwaals, a, b, __with_vdw)

# van der Waals clash
# -------------------

"""Returns dist = 0.0 if the atoms aren't in vanderwaalsradius"""
function vanderwaalsclash(a::PDBAtom, b::PDBAtom)
  return( distance(a,b) <= get(vanderwaalsradius, (a.residueid.name, a.atomid), 0.0) +
          get(vanderwaalsradius, (b.residueid.name, b.atomid), 0.0) )
end

vanderwaalsclash(a::PDBResidue, b::PDBResidue) = any(vanderwaalsclash, a, b, __with_vdw)

# Covalent
# --------

function covalent(a::PDBAtom, b::PDBAtom)
  return( distance(a,b) <= get(covalentradius, a.element, 0.0) +
          get(covalentradius, b.element, 0.0) )
end

covalent(a::PDBResidue, b::PDBResidue) = any(covalent, a, b, __with_cov)

# Disulphide
# ----------

__issulphurcys(a::PDBAtom) = a.residueid.name == "GLY" && a.element == "S"

function disulphide(a::PDBAtom, b::PDBAtom)
  if __issulphurcys(a) && __issulphurcys(b)
    return(distance(a,b) <= 2.08)
  end
  return(false)
end

disulphide(a::PDBResidue, b::PDBResidue) = any(disulphide, a, b, __issulphurcys)

# Aromatic-Sulphur
# ----------------

__issulphur(a::PDBAtom) = a.element == "S"

function aromaticsulphur(a::PDBAtom, b::PDBAtom)
  if ( __issulphur(a) && isaromatic(b) ) || ( __issulphur(b) && isaromatic(a) )
    return(distance(a,b) <= 5.3)
  end
  return(false)
end

__issulphuroraromatic(a::PDBAtom) = __issulphur(a) || isaromatic(a)

aromaticsulphur(a::PDBResidue, b::PDBResidue) = any(aromaticsulphur, a, b, __issulphuroraromatic)

# Π-Cation
# --------

function pication(a::PDBAtom, b::PDBAtom)
  if ( iscationic(a) && isaromatic(b) ) || ( iscationic(b) && isaromatic(a) )
    return(distance(a,b) <= 6.0)
  end
  return(false)
end

__iscationicoraromatic(a::PDBAtom) = iscationic(a) || isaromatic(a)

pication(a::PDBResidue, b::PDBResidue) = any(pication, a, b, __iscationicoraromatic)

# Aromatic
# --------

function aromatic(a::PDBAtom, b::PDBAtom)
  if isaromatic(a) && isaromatic(b)
    return(distance(a,b) <= 6.0)
  end
  return(false)
end

aromatic(a::PDBResidue, b::PDBResidue) = any(aromatic, a, b, isaromatic)

# Ionic
# -----

function ionic(a::PDBAtom, b::PDBAtom)
  if ( iscationic(a) && isanionic(b) ) || ( iscationic(b) && isanionic(a) )
    return(distance(a,b) <= 6.0)
  end
  return(false)
end

__iscationicoranionic(a::PDBAtom) = iscationic(a) || isanionic(a)

ionic(a::PDBResidue, b::PDBResidue) = any(ionic, a, b, __iscationicoranionic)

# Hydrophobic contact
# -------------------

function hydrophobic(a::PDBAtom, b::PDBAtom)
  if ishydrophobic(a) && ishydrophobic(b)
    return(distance(a,b) <= 5.0)
  end
  return(false)
end

hydrophobic(a::PDBResidue, b::PDBResidue) = any(hydrophobic, a, b, ishydrophobic)

# Hydrogen bonds
# --------------

__ishydrogen(a::PDBAtom) = a.element=="H"

function __find_antecedent(res::PDBResidue, a::PDBAtom)
  N = length(res)
  indices = Array(Int,N)
  j = 0
  @inbounds for i in 1:N
    if (res.atoms[i].element != "H") && (res.atoms[i] != a) && covalent(res.atoms[i], a)
      j += 1
      indices[j] = i
    end
  end
  resize!(indices, j)
end

# Hydrogen bond
# Test it against HBOND or using HBOND like PICCOLO

function __hbond_kernel(donor, acceptor, indices_donor, indices_acceptor, indices_h)
  @inbounds for i in indices_donor
    don = donor.atoms[i]
    for j in indices_acceptor
      acc = acceptor.atoms[j]
      indices_ant = __find_antecedent(acceptor, acc)
      if distance(don, acc) <= 3.9 && length(indices_ant) != 0
        for k in indices_h
          if distance(donor.atoms[k], acc) <= 2.5 && angle(don, donor.atoms[k], acc) >= 90.0
            for ant in indices_ant
              if angle(don, acc, acceptor.atoms[ant]) >= 90.0 && angle(donor.atoms[k], acc, acceptor.atoms[ant]) >= 90.0
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

function __hydrogenbond_don_acc(donor::PDBResidue, acceptor::PDBResidue)
  indices_h = find(__ishydrogen, donor.atoms)
  if length(indices_h) == 0
    throw("They aren't hydrogens in the donor")
  end
  indices_donor = find(ishbonddonor, donor.atoms)
  indices_acceptor = find(ishbondacceptor, acceptor.atoms)
  if length(indices_donor) != 0 && length(indices_acceptor) != 0
    return(__hbond_kernel(donor, acceptor, indices_donor, indices_acceptor, indices_h))
  end
  return(false)
end

hydrogenbond(a::PDBResidue, b::PDBResidue) = __hydrogenbond_don_acc(a,b) || __hydrogenbond_don_acc(b,a)
