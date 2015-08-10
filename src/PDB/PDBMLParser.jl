using LightXML

function _get_text(elem, name)
	sub = find_element(elem, name)
	if sub != nothing
		return(content(sub))
	end
	throw("There is not $name for $elem")
end

function _get_atom_iterator(pdbml::ASCIIString)
	pdbfile = parse_file(pdbml)
	pdbroot = root(pdbfile)
	child_elements(get_elements_by_tagname(pdbroot, "atom_siteCategory")[1])
end

function getpdbmlatoms(pdbml::ASCIIString; chain::ASCIIString = "all",
	model::ASCIIString = "all", group::ASCIIString = "all", atomname::ASCIIString="all", onlyheavy::Bool=false)
	atom_list = Array(PDBAtom,0)
	atoms = _get_atom_iterator(pdbml)
	for atom in atoms

   	atom_group = _get_text(atom, "group_PDB")
		atom_model = _get_text(atom, "pdbx_PDB_model_num")
		atom_chain = _get_text(atom, "label_asym_id")
   	atom_name = _get_text(atom, "label_atom_id")
   	element = _get_text(atom, "type_symbol")

   	if  (group=="all" || group==atom_group) && (chain=="all" || chain==atom_chain) &&
   	    (model=="all" || model==atom_model) && (atomname=="all" || atomname==atom_name) && (!onlyheavy || element!="H")

			number = _get_text(atom, "label_seq_id")
			name = _get_text(atom, "label_comp_id")
			x = float(_get_text(atom, "Cartn_x"))
			y = float(_get_text(atom, "Cartn_y"))
			z = float(_get_text(atom, "Cartn_z"))
			occupancy = float(_get_text(atom, "occupancy"))
			B = _get_text(atom, "B_iso_or_equiv")

			push!(atom_list, PDBAtom(PDBResidueIdentifier(number, name, atom_group, atom_model, atom_chain),
                               Coordinates(x,y,z), atom_name, element, occupancy, B))
		end
	end
	atom_list
end

function _add_atom!(res::PDBResidue, atom::PDBAtom)
	if res.id == atom.residueid
		push!(res.atoms, atom)
		return(res)
	else
		throw("It isn't the same residue: $(res.id) != $(atom.residueid) ")
	end
end

function getresidues(atoms::Vector{PDBAtom})
  residues = Dict{PDBResidueIdentifier, PDBResidue}()
  for atom in atoms
    if atom.residueid in keys(residues)
      residues[atom.residueid] = _add_atom!(residues[atom.residueid], atom)
    else
      residues[atom.residueid] = PDBResidue(atom.residueid, [atom])
    end
  end
  sizehint!(residues, length(residues))
end

getresidues(pdbml::ASCIIString; chain::ASCIIString = "all",	model::ASCIIString = "all",
            group::ASCIIString = "all", onlyheavy::Bool=false) = getresidues(getatoms(pdbml, chain=chain, model=model, group=group, onlyheavy=onlyheavy))

# Download PDB
# ============

function _inputnameforgzip(outfile)
  len = length(outfile)
  if len > 3 && outfile[len-2:len] == ".gz"
    return(outfile)
  end
  string(outfile, ".gz")
end

function downloadpdb(pdbcode::ASCIIString; format::ASCIIString="xml", outfile::ASCIIString="default")
  if length(pdbcode)== 4
    filename = string(uppercase(pdbcode), ".", lowercase(format),".gz")
    outfile = outfile == "default" ? filename : _inputnameforgzip(outfile)
    download(string("http://www.rcsb.org/pdb/files/", filename), outfile)
    run(`gzip -d $outfile`)
  else
    throw(string(pdbcode, " is not a correct PDB code"))
  end
end
