using LightXML

function __get_text(elem, name)
	sub = find_element(elem, name)
	if sub != nothing
		return(content(sub))
	end
	throw("There is not $name for $elem")
end

function getpdbchain(filename::ASCIIString; inchain::ASCIIString = "A",
	inmodel::ASCIIString = "1", ingroup::ASCIIString = "ATOM")
	pdbchain = Dict{ASCIIString,PDBResidue}()
	pdbfile = parse_file(filename);
	pdbroot = root(pdbfile);
	for atom in child_elements(get_elements_by_tagname(pdbroot, "atom_siteCategory")[1])
		group = __get_text(atom, "group_PDB")
		model = __get_text(atom, "pdbx_PDB_model_num")
		chain = __get_text(atom, "auth_asym_id")
		if ingroup==group && inchain==chain && inmodel==model
			number = __get_text(atom, "auth_seq_id")
			name = __get_text(atom, "auth_comp_id")
			coordinates = ( float(__get_text(atom, "Cartn_x")), float(__get_text(atom, "Cartn_y")), float(__get_text(atom, "Cartn_z")) )
			atom_name = __get_text(atom, "auth_atom_id")
			element = __get_text(atom, "type_symbol")
			occupancy = float(__get_text(atom, "occupancy"))
			B = __get_text(atom, "B_iso_or_equiv")
			if number in keys(pdbchain)
				pdbchain[number] = __add_atom!(pdbchain[number], number, name, group, atom_name, coordinates, element, occupancy, B)
			else
				pdbchain[number] = pdbresidue(number, name, group, atom_name, coordinates, element, occupancy, B)
			end
		end
	end
	pdbchain
end