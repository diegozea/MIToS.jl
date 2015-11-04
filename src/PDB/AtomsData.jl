"""Covalent radius in Å of each element from the Additional file 1 of PICCOLO [1].
Hydrogen was updated using the value on Table 2 from Cordero et. al. [2].
1. Bickerton, G. R., Higueruelo, A. P., & Blundell, T. L. (2011).
Comprehensive, atomic-level characterization of structurally
characterized protein-protein interactions: the PICCOLO database.
BMC bioinformatics, 12(1), 313.
2. Cordero, B., Gómez, V., Platero-Prats, A. E., Revés, M.,
Echeverría, J., Cremades, E., ... & Alvarez, S. (2008).
Covalent radii revisited. Dalton Transactions, (21), 2832-2838."""
const covalentradius = Dict{ASCIIString,Float64}("C" => 0.77,
                                                 "N" => 0.70,
                                                 "O" => 0.66,
                                                 "S" => 1.04,
                                                 "H" => 0.31 )

const _3_letter_aa = ASCIIString[ "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL" ]

"""van der Waals radius in Å from the Additional file 1 of
Bickerton, G. R., Higueruelo, A. P., & Blundell, T. L. (2011).
Comprehensive, atomic-level characterization of structurally characterized protein-protein interactions: the PICCOLO database.
BMC bioinformatics, 12(1), 313."""
const vanderwaalsradius = Dict{Tuple{ASCIIString,ASCIIString},Float64}(("ALA","C") => 1.61,
("ALA","CA") => 1.88,
("ALA","CB") => 1.88,
("ALA","N") => 1.64,
("ALA","O") => 1.42,
("ARG","C") => 1.61,
("ARG","CA") => 1.88,
("ARG","CB") => 1.88,
("ARG","CD") => 1.88,
("ARG","CG") => 1.88,
("ARG","CZ") => 1.61,
("ARG","N") => 1.64,
("ARG","NE") => 1.64,
("ARG","NH1") => 1.64,
("ARG","NH2") => 1.64,
("ARG","O") => 1.42,
("ASN","C") => 1.61,
("ASN","CA") => 1.88,
("ASN","CB") => 1.88,
("ASN","CG") => 1.61,
("ASN","N") => 1.64,
("ASN","ND2") => 1.64,
("ASN","O") => 1.42,
("ASN","OD1") => 1.42,
("ASP","C") => 1.61,
("ASP","CA") => 1.88,
("ASP","CB") => 1.88,
("ASP","CG") => 1.61,
("ASP","N") => 1.64,
("ASP","O") => 1.42,
("ASP","OD1") => 1.42,
("ASP","OD2") => 1.42,
("CYS","C") => 1.61,
("CYS","CA") => 1.88,
("CYS","CB") => 1.88,
("CYS","N") => 1.64,
("CYS","O") => 1.42,
("CYS","SG") => 1.77,
("GLN","C") => 1.61,
("GLN","CA") => 1.88,
("GLN","CB") => 1.88,
("GLN","CD") => 1.61,
("GLN","CG") => 1.88,
("GLN","N") => 1.64,
("GLN","NE2") => 1.64,
("GLN","O") => 1.42,
("GLN","OE1") => 1.42,
("GLU","C") => 1.61,
("GLU","CA") => 1.88,
("GLU","CB") => 1.88,
("GLU","CD") => 1.61,
("GLU","CG") => 1.88,
("GLU","N") => 1.64,
("GLU","O") => 1.42,
("GLU","OE1") => 1.42,
("GLU","OE2") => 1.42,
("GLY","C") => 1.61,
("GLY","CA") => 1.88,
("GLY","N") => 1.64,
("GLY","O") => 1.42,
("HIS","C") => 1.61,
("HIS","CA") => 1.88,
("HIS","CB") => 1.88,
("HIS","CD2") => 1.76,
("HIS","CE1") => 1.76,
("HIS","CG") => 1.61,
("HIS","N") => 1.64,
("HIS","ND1") => 1.64,
("HIS","NE2") => 1.64,
("HIS","O") => 1.42,
("ILE","C") => 1.61,
("ILE","CA") => 1.88,
("ILE","CB") => 1.88,
("ILE","CD1") => 1.88,
("ILE","CG1") => 1.88,
("ILE","CG2") => 1.88,
("ILE","N") => 1.64,
("ILE","O") => 1.42,
("LEU","C") => 1.61,
("LEU","CA") => 1.88,
("LEU","CB") => 1.88,
("LEU","CD1") => 1.88,
("LEU","CD2") => 1.88,
("LEU","CG") => 1.88,
("LEU","N") => 1.64,
("LEU","O") => 1.42,
("LYS","C") => 1.61,
("LYS","CA") => 1.88,
("LYS","CB") => 1.88,
("LYS","CD") => 1.88,
("LYS","CE") => 1.88,
("LYS","CG") => 1.88,
("LYS","N") => 1.64,
("LYS","NZ") => 1.64,
("LYS","O") => 1.42,
("MET","C") => 1.61,
("MET","CA") => 1.88,
("MET","CB") => 1.88,
("MET","CE") => 1.88,
("MET","CG") => 1.88,
("MET","N") => 1.64,
("MET","O") => 1.42,
("MET","SD") => 1.77,
("PHE","C") => 1.61,
("PHE","CA") => 1.88,
("PHE","CB") => 1.88,
("PHE","CD1") => 1.76,
("PHE","CD2") => 1.76,
("PHE","CE1") => 1.76,
("PHE","CE2") => 1.76,
("PHE","CG") => 1.61,
("PHE","CZ") => 1.76,
("PHE","N") => 1.64,
("PHE","O") => 1.42,
("PRO","C") => 1.61,
("PRO","CA") => 1.88,
("PRO","CB") => 1.88,
("PRO","CD") => 1.88,
("PRO","CG") => 1.88,
("PRO","N") => 1.64,
("PRO","O") => 1.42,
("SER","C") => 1.61,
("SER","CA") => 1.88,
("SER","CB") => 1.88,
("SER","N") => 1.64,
("SER","O") => 1.42,
("SER","OG") => 1.46,
("THR","C") => 1.61,
("THR","CA") => 1.88,
("THR","CB") => 1.88,
("THR","CG2") => 1.88,
("THR","N") => 1.64,
("THR","O") => 1.42,
("THR","OG1") => 1.46,
("TRP","C") => 1.61,
("TRP","CA") => 1.88,
("TRP","CB") => 1.88,
("TRP","CD1") => 1.76,
("TRP","CD2") => 1.61,
("TRP","CE2") => 1.61,
("TRP","CE3") => 1.76,
("TRP","CG") => 1.61,
("TRP","CH2") => 1.76,
("TRP","CZ2") => 1.76,
("TRP","CZ3") => 1.76,
("TRP","N") => 1.64,
("TRP","NE1") => 1.64,
("TRP","O") => 1.42,
("TYR","C") => 1.61,
("TYR","CA") => 1.88,
("TYR","CB") => 1.88,
("TYR","CD1") => 1.76,
("TYR","CD2") => 1.76,
("TYR","CE1") => 1.76,
("TYR","CE2") => 1.76,
("TYR","CG") => 1.61,
("TYR","CZ") => 1.61,
("TYR","N") => 1.64,
("TYR","O") => 1.42,
("TYR","OH") => 1.46,
("VAL","C") => 1.61,
("VAL","CA") => 1.88,
("VAL","CB") => 1.88,
("VAL","CG1") => 1.88,
("VAL","CG2") => 1.88,
("VAL","N") => 1.64,
("VAL","O") => 1.42 )

function _add_CTER_O!(dict)
  for aa in _3_letter_aa
    push!(dict, (aa, "OXT"))
    push!(dict, (aa, "OT2"))
    push!(dict, (aa, "OT1"))
  end
  dict
end

function _add_CTER_O!(dict, value)
  for aa in _3_letter_aa
    push!(dict, (aa, "OXT") => value)
    push!(dict, (aa, "OT2") => value)
    push!(dict, (aa, "OT1") => value)
  end
  dict
end

_add_CTER_O!(vanderwaalsradius, 1.42) # Using 1.42 because OXT is the terminal oxigen I assume same vdw radii than O

const _hydrophobic = Set{Tuple{ASCIIString, ASCIIString}}( [ ("ALA","CB"),
("ARG","CB"),
("ARG","CG"),
("ASN","CB"),
("ASP","CB"),
("CYS","CB"),
("GLN","CB"),
("GLN","CG"),
("GLU","CB"),
("GLU","CG"),
("HIS","CB"),
("ILE","CB"),
("ILE","CD1"),
("ILE","CG1"),
("ILE","CG2"),
("LEU","CB"),
("LEU","CD1"),
("LEU","CD2"),
("LEU","CG"),
("LYS","CB"),
("LYS","CD"),
("LYS","CG"),
("MET","CB"),
("MET","CE"),
("MET","CG"),
("MET","SD"),
("PHE","CB"),
("PHE","CD1"),
("PHE","CD2"),
("PHE","CE1"),
("PHE","CE2"),
("PHE","CG"),
("PHE","CZ"),
("PRO","CB"),
("PRO","CG"),
("THR","CG2"),
("TRP","CB"),
("TRP","CD2"),
("TRP","CE3"),
("TRP","CG"),
("TRP","CH2"),
("TRP","CZ2"),
("TRP","CZ3"),
("TYR","CB"),
("TYR","CD1"),
("TYR","CD2"),
("TYR","CE1"),
("TYR","CE2"),
("TYR","CG"),
("VAL","CB"),
("VAL","CG1"),
("VAL","CG2") ] )

const _aromatic_res = Set{ASCIIString}( [ "HIS", "PHE", "TRP", "TYR" ] )

const _aromatic = Set{Tuple{ASCIIString, ASCIIString}}( [ ("HIS","CD2"),
("HIS","CE1"),
("HIS","CG"),
("HIS","ND1"),
("HIS","NE2"),
("PHE","CD1"),
("PHE","CD2"),
("PHE","CE1"),
("PHE","CE2"),
("PHE","CG"),
("PHE","CZ"),
("TRP","CD1"),
("TRP","CD2"),
("TRP","CE2"),
("TRP","CE3"),
("TRP","CG"),
("TRP","CH2"),
("TRP","CZ2"),
("TRP","CZ3"),
("TRP","NE1"),
("TYR","CD1"),
("TYR","CD2"),
("TYR","CE1"),
("TYR","CE2"),
("TYR","CG"),
("TYR","CZ") ] )

const _cationic = Set{Tuple{ASCIIString, ASCIIString}}( [ ("ARG","CZ"),
("ARG","NE"),
("ARG","NH1"),
("ARG","NH2"),
("HIS","CD2"),
("HIS","CE1"),
("HIS","CG"),
("HIS","ND1"),
("HIS","NE2"),
("LYS","NZ") ] )

const _anionic = Set{Tuple{ASCIIString, ASCIIString}}( [ ("ASP","CG"),
("ASP","OD1"),
("ASP","OD2"),
("GLU","CD"),
("GLU","OE1"),
("GLU","OE2") ] )

_add_CTER_O!(_anionic)

"""Keys come from Table 1 of Bickerton et. al. 2011,
The hydrogen names of the donor comes from: http://biomachina.org/courses/modeling/download/topallh22x.pro
Synonyms come from: http://www.bmrb.wisc.edu/ref_info/atom_nom.tbl"""
const _hbond_donor = Dict{Tuple{ASCIIString, ASCIIString}, Vector{ASCIIString}}( ("ALA","N") => ["HN", "H", "HN1","H1","1H", "HN2","H2","2H", "HN3","H3","3H", "HT1","HT2","HT3"],
("ARG","N") => ["HN", "H", "HN1","H1","1H", "HN2","H2","2H", "HN3","H3","3H", "HT1","HT2","HT3"],
("ARG","NE") => ["HE", "HNE"],
("ARG","NH1") => ["HH11","HH12", "1HH1","2HH1"],
("ARG","NH2") => ["HH22","HH21", "2HH1", "1HH2"],
("ASN","N") => ["HN", "H", "HN1","H1","1H", "HN2","H2","2H", "HN3","H3","3H", "HT1","HT2","HT3"],
("ASN","ND2") => ["HD21","HD22", "HN21", "HN22", "1HD2","2HD2"],
("ASP","N") => ["HN", "H", "HN1","H1","1H", "HN2","H2","2H", "HN3","H3","3H", "HT1","HT2","HT3"],
("CYS","N") => ["HN", "H", "HN1","H1","1H", "HN2","H2","2H", "HN3","H3","3H", "HT1","HT2","HT3"],
("CYS","SG") => ["HG1", "HG", "HSG"],
("GLN","N") => ["HN", "H", "HN1","H1","1H", "HN2","H2","2H", "HN3","H3","3H", "HT1","HT2","HT3"],
("GLN","NE2") => ["HE21","HE22", "1HE2","2HE2", "HN21","HN22"],
("GLU","N") => ["HN", "H", "HN1","H1","1H", "HN2","H2","2H", "HN3","H3","3H", "HT1","HT2","HT3"],
("GLY","N") => ["HN", "H", "HN1","H1","1H", "HN2","H2","2H", "HN3","H3","3H", "HT1","HT2","HT3"],
("HIS","N") => ["HN", "H", "HN1","H1","1H", "HN2","H2","2H", "HN3","H3","3H", "HT1","HT2","HT3"],
("HIS","ND1") => ["HD1"],
("HIS","NE2") => [ "HE2","HNE2"],
("ILE","N") => ["HN", "H", "HN1","H1","1H", "HN2","H2","2H", "HN3","H3","3H", "HT1","HT2","HT3"],
("LEU","N") => ["HN", "H", "HN1","H1","1H", "HN2","H2","2H", "HN3","H3","3H", "HT1","HT2","HT3"],
("LYS","N") => ["HN", "H", "HN1","H1","1H", "HN2","H2","2H", "HN3","H3","3H", "HT1","HT2","HT3"],
("LYS","NZ") => ["HZ1","HZ2","HZ3", "1HZ","2HZ","3HZ", "HNZ1","HNZ2","HNZ3"],
("MET","N") => ["HN", "H", "HN1","H1","1H", "HN2","H2","2H", "HN3","H3","3H", "HT1","HT2","HT3"],
("PHE","N") => ["HN", "H", "HN1","H1","1H", "HN2","H2","2H", "HN3","H3","3H", "HT1","HT2","HT3"],
("SER","N") => ["HN", "H", "HN1","H1","1H", "HN2","H2","2H", "HN3","H3","3H", "HT1","HT2","HT3"],
("SER","OG") => ["HG1", "HG", "HOG"],
("THR","N") => ["HN", "H", "HN1","H1","1H", "HN2","H2","2H", "HN3","H3","3H", "HT1","HT2","HT3"],
("THR","OG1") => ["HG1", "HOG1"],
("TRP","N") => ["HN", "H", "HN1","H1","1H", "HN2","H2","2H", "HN3","H3","3H", "HT1","HT2","HT3"],
("TRP","NE1") => ["HE1", "HNE1"],
("TYR","N") => ["HN", "H", "HN1","H1","1H", "HN2","H2","2H", "HN3","H3","3H", "HT1","HT2","HT3"],
("TYR","OH") => ["HH", "HOH"],
("VAL","N") => ["HN", "H", "HN1","H1","1H", "HN2","H2","2H", "HN3","H3","3H", "HT1","HT2","HT3"],
("PRO","N") => ["HN1","H1","1H", "HN2","H2","2H", "HT1","HT2"] )
# Proline N-Terminal (RESIDUE PROP)


"""Keys come from Table 1 of Bickerton et. al. 2011,
Antecedents come from come from: http://biomachina.org/courses/modeling/download/topallh22x.pro
Synonyms come from: http://www.bmrb.wisc.edu/ref_info/atom_nom.tbl"""
const _hbond_acceptor = Dict{Tuple{ASCIIString, ASCIIString}, Vector{ASCIIString}}( ("ALA","O") => ["C"], ("ALA","OT1") => ["C"], ("ALA","OXT") => ["C"], ("ALA","OT2") => ["C"],
("ARG","O") => ["C"], ("ARG","OT1") => ["C"], ("ARG","OXT") => ["C"], ("ARG","OT2") => ["C"],
("ASN","O") => ["C"], ("ASN","OT1") => ["C"], ("ASN","OXT") => ["C"], ("ASN","OT2") => ["C"],
("ASN","OD1") => ["CG"],
("ASP","O") => ["C"], ("ASP","OT1") => ["C"], ("ASP","OXT") => ["C"], ("ASP","OT2") => ["C"],
("ASP","OD1") => ["CG"],
("ASP","OD2") => ["CG"],
("CYS","O") => ["C"], ("CYS","OT1") => ["C"], ("CYS","OXT") => ["C"], ("CYS","OT2") => ["C"],
("CYS","SG") => ["CB"],
("GLN","O") => ["C"], ("GLN","OT1") => ["C"], ("GLN","OXT") => ["C"], ("GLN","OT2") => ["C"],
("GLN","OE1") => ["CD"],
("GLU","O") => ["C"], ("GLU","OT1") => ["C"], ("GLU","OXT") => ["C"], ("GLU","OT2") => ["C"],
("GLU","OE1") => ["CD"],
("GLU","OE2") => ["CD"],
("GLY","O") => ["C"], ("GLY","OT1") => ["C"], ("GLY","OXT") => ["C"], ("GLY","OT2") => ["C"],
("HIS","ND1") => ["CG", "CE1"],
("HIS","NE2") => ["CD2", "CE1"],
("HIS","O") => ["C"], ("HIS","OT1") => ["C"], ("HIS","OXT") => ["C"], ("HIS","OT2") => ["C"],
("ILE","O") => ["C"], ("ILE","OT1") => ["C"], ("ILE","OXT") => ["C"], ("ILE","OT2") => ["C"],
("LEU","O") => ["C"], ("LEU","OT1") => ["C"], ("LEU","OXT") => ["C"], ("LEU","OT2") => ["C"],
("LYS","O") => ["C"], ("LYS","OT1") => ["C"], ("LYS","OXT") => ["C"], ("LYS","OT2") => ["C"],
("MET","O") => ["C"], ("MET","OT1") => ["C"], ("MET","OXT") => ["C"], ("MET","OT2") => ["C"],
("MET","SD") => ["CG", "CE"],
("PHE","O") => ["C"], ("PHE","OT1") => ["C"], ("PHE","OXT") => ["C"], ("PHE","OT2") => ["C"],
("PRO","O") => ["C"], ("PRO","OT1") => ["C"], ("PRO","OXT") => ["C"], ("PRO","OT2") => ["C"],
("SER","O") => ["C"], ("SER","OT1") => ["C"], ("SER","OXT") => ["C"], ("SER","OT2") => ["C"],
("SER","OG") => ["CB"],
("THR","O") => ["C"], ("THR","OT1") => ["C"], ("THR","OXT") => ["C"], ("THR","OT2") => ["C"],
("THR","OG1") => ["CB"],
("TRP","O") => ["C"], ("TRP","OT1") => ["C"], ("TRP","OXT") => ["C"], ("TRP","OT2") => ["C"],
("TYR","O") => ["C"], ("TYR","OT1") => ["C"], ("TYR","OXT") => ["C"], ("TYR","OT2") => ["C"],
("VAL","O") => ["C"], ("VAL","OT1") => ["C"], ("VAL","OXT") => ["C"], ("VAL","OT2") => ["C"] )

function _generate_dict!(dict, input_dict)
  for (res, atom) in keys(input_dict)
    if haskey(dict, res)
      push!(dict[res], atom)
    else
      dict[res] = Set{ASCIIString}(ASCIIString[ atom ])
    end
  end
  dict
end

function _generate_dict!(dict, input_set::Set{Tuple{ASCIIString,ASCIIString}})
  for (res, atom) in input_set
    if haskey(dict, res)
      push!(dict[res], atom)
    else
      dict[res] = Set{ASCIIString}(ASCIIString[ atom ])
    end
  end
  dict
end

function _generate_interaction_keys(vdw, hyd, aro, cat, ani)
  dict = Dict{ASCIIString, Set{ASCIIString}}()
  _generate_dict!(dict, vdw)
  _generate_dict!(dict, hyd)
  _generate_dict!(dict, aro)
  _generate_dict!(dict, cat)
  _generate_dict!(dict, ani)
  dict
end

const _interaction_keys = _generate_interaction_keys(vanderwaalsradius, _hydrophobic, _aromatic, _cationic, _anionic)

_generate_atoms_set(res::PDBResidue) = ASCIIString[ atom.atom for atom in res.atoms[findheavy(res)] ]

function check_atoms_for_interactions(res::PDBResidue)
  atoms = _generate_atoms_set(res)
  if haskey(_interaction_keys, res.id.name)
    used = _interaction_keys[res.id.name]
  else
    warn(string("RESIDUE: ", res.id.name, " is unknown for MIToS.PDB (AtomsData.jl)"))
    return(false)
  end
  for atom in atoms
    if !( atom in used )
      warn(string("RESIDUE ", res.id.name," ATOM ", atom, " is unknown for MIToS.PDB (AtomsData.jl)"))
      return( false )
    end
  end
  true
end
