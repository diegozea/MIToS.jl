"""Group can be ATOM or HETATM"""
function getpdbatoms(pdb::ASCIIString; chain::ASCIIString = "all",
                     model::ASCIIString = "all", group::ASCIIString = "all", atomname::ASCIIString="all", onlyheavy::Bool=false)
  residue_dict = OrderedDict{PDBResidueIdentifier, Vector{PDBAtom}}()
  fh = open(pdb, "r")
  atom_model = 0
  for line in eachline(fh)
    line_id = replace(line[1:6], ' ', "")
    if line_id == "MODEL"
      atom_model += 1
    end
    if (group=="all" && (line_id=="ATOM" || line_id=="HETATM")) || (line_id == group)
      atom_chain = string(line[22])
      atom_name = replace(line[13:16],' ',"")
      element = replace(line[77:78],' ',"")
      if  (chain=="all" || chain==atom_chain) && ((model==0 || model=="all") || model==string(atom_model+1)) &&
          (atomname=="all" || atomname==atom_name) && (!onlyheavy || element!="H")

        # 23 - 26        Integer         Residue sequence number.
        # 27             AChar           Code for insertion of residues.
        PDB_number = replace(line[23:27],' ',"")

        name = replace(line[18:20],' ',"")
        x = float(replace(line[31:38],' ',""))
        y = float(replace(line[39:46],' ',""))
        z = float(replace(line[47:54],' ',""))
        occupancy = float(replace(line[55:60],' ',""))
        B = replace(line[61:66],' ',"")

        mdl = atom_model == 0 ? "1" : string(atom_model)

        residue_id = PDBResidueIdentifier("", PDB_number, name, line_id, mdl, atom_chain)
        atom_data  = PDBAtom(Coordinates(x,y,z), atom_name, element, occupancy, B)

        value = get!(residue_dict, residue_id, PDBAtom[])
        push!(value, atom_data)

      end
    end
  end
  close(fh)
  _generate_residues(residue_dict)
end
