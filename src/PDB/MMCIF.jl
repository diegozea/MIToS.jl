struct MMCIFFile <: FileFormat end

function _parse_mmcif_to_pdbresidues(mmcif_dict::MMCIFDict, label::Bool)
    # Choose the correct prefix based on the label argument
    prefix = label ? "_atom_site.label" : "_atom_site.auth"
    chain_attr = string(prefix, "_asym_id")
    comp_id_attr = string(prefix, "_comp_id")
    atom_id_attr = string(prefix, "_atom_id")

    # Extract relevant entries from the MMCIFDict
    auth_asym_ids = mmcif_dict[chain_attr] # Chain
    auth_seq_ids = mmcif_dict["_atom_site.auth_seq_id"] # Residue number
    label_seq_ids = mmcif_dict["_atom_site.label_seq_id"] # Residue name (PDBe)
    auth_comp_ids = mmcif_dict[comp_id_attr] # Residue name (PDB)
    atom_names = mmcif_dict[atom_id_attr] # Atom name, e.g. "CA"
    cartn_x = mmcif_dict["_atom_site.Cartn_x"] # x
    cartn_y = mmcif_dict["_atom_site.Cartn_y"] # y
    cartn_z = mmcif_dict["_atom_site.Cartn_z"] # z
    occupancies = mmcif_dict["_atom_site.occupancy"] # Occupancy
    bfactors = mmcif_dict["_atom_site.B_iso_or_equiv"] # B-factor
    elements = mmcif_dict["_atom_site.type_symbol"] # Element, e.g. "C"
    group_pdb = mmcif_dict["_atom_site.group_PDB"]  # Group type, "ATOM" or "HETATM"
    pdb_model = mmcif_dict["_atom_site.pdbx_PDB_model_num"]  # Model number

    ins_codes = mmcif_dict["_atom_site.pdbx_PDB_ins_code"] # Insertion codes

    residues = PDBResidue[]
    current_residue_id = ""
    current_residue = PDBResidue(PDBResidueIdentifier("", "", "", "", "", ""), PDBAtom[])

    for i = 1:length(auth_seq_ids)
        pdb_number = string(auth_seq_ids[i], replace(ins_codes[i], "?" => ""))
        pdbe_number = replace(label_seq_ids[i], "." => "")
        residue_id = PDBResidueIdentifier(
            pdbe_number,
            pdb_number,
            auth_comp_ids[i],
            group_pdb[i],
            pdb_model[i],
            auth_asym_ids[i],
        )

        if current_residue_id != string(residue_id)
            if !isempty(current_residue.atoms)
                push!(residues, current_residue)
            end
            current_residue = PDBResidue(residue_id, Vector{PDBAtom}())
            current_residue_id = string(residue_id)
        end

        atom = PDBAtom(
            Coordinates(
                parse(Float64, cartn_x[i]),
                parse(Float64, cartn_y[i]),
                parse(Float64, cartn_z[i]),
            ),
            atom_names[i],
            elements[i],
            parse(Float64, occupancies[i]),
            bfactors[i],
        )

        push!(current_residue.atoms, atom)
    end

    if !isempty(current_residue.atoms)
        push!(residues, current_residue)
    end

    residues
end

"""
`parse_file(io, ::Type{MMCIFFile}; chain=All, model=All, group=All, atomname=All, onlyheavy=false, label=true, occupancyfilter=false)`

Parse an mmCIF file and returns a list of `PDBResidue`s. Setting `chain`, `model`, `group`,
`atomname` and `onlyheavy` values can be used to select a subset of residues. Group can be
`"ATOM"` or `"HETATM"`. If those keyword arguments are not set, all residues are returned.
If the keyword argument `label` (default: `true`) is `false`, the **auth_** attributes will be used instead
of the **label_** attributes for `chain`, `atom`, and residue `name` fields. The **auth_**
attributes are alternatives provided by an author in order to match the
identification/values used in the publication that describes the structure. If the
keyword argument `occupancyfilter` (default: `false`) is `true`, only the atoms
with the best occupancy are returned.
"""
function Utils.parse_file(
    io::Union{IO,String},
    ::Type{MMCIFFile};
    chain::Union{String,Type{All}} = All,
    model::Union{String,Type{All}} = All,
    group::Union{String,Type{All}} = All,
    atomname::Union{String,Type{All}} = All,
    onlyheavy::Bool = false,
    label::Bool = true,
    occupancyfilter::Bool = false,
)
    mmcif_dict = MMCIFDict(io)

    residues = select_residues(
        _parse_mmcif_to_pdbresidues(mmcif_dict, label),
        model = model,
        chain = chain,
        group = group,
    )

    for res in residues
        filter!(a -> _is(a.atom, atomname), res.atoms)
        if occupancyfilter
            res.atoms = bestoccupancy(res.atoms)
        end
        if onlyheavy
            filter!(a -> a.element != "H", res.atoms)
        end
    end

    filter!(res -> !isempty(res.atoms), residues)
end

function _resnumber(number)
    num = replace(number, r"[A-Za-z]" => "")
    isempty(num) ? "." : num
end

function _inscode(res::PDBResidue)
    m = match(r"[A-Za-z]$", res.id.number)
    return m === nothing ? "?" : m.match
end

function _pdbresidues_to_mmcifdict(
    residues::Vector{PDBResidue};
    label::Bool = false,
    molecular_structures::Bool = false,
)
    # Initialize MMCIFDict with the necessary fields
    mmcif_dict = MMCIFDict()

    # Initialize fields as empty arrays
    if molecular_structures || !label
        mmcif_dict["_atom_site.auth_asym_id"] = String[]
        mmcif_dict["_atom_site.auth_comp_id"] = String[]
        mmcif_dict["_atom_site.auth_atom_id"] = String[]
    end
    if label
        mmcif_dict["_atom_site.label_asym_id"] = String[]
        mmcif_dict["_atom_site.label_comp_id"] = String[]
        mmcif_dict["_atom_site.label_atom_id"] = String[]
    end
    mmcif_dict["_atom_site.id"] = String[]
    mmcif_dict["_atom_site.auth_seq_id"] = String[]
    mmcif_dict["_atom_site.label_seq_id"] = String[]
    mmcif_dict["_atom_site.Cartn_x"] = String[]
    mmcif_dict["_atom_site.Cartn_y"] = String[]
    mmcif_dict["_atom_site.Cartn_z"] = String[]
    mmcif_dict["_atom_site.occupancy"] = String[]
    mmcif_dict["_atom_site.B_iso_or_equiv"] = String[]
    mmcif_dict["_atom_site.type_symbol"] = String[]
    mmcif_dict["_atom_site.group_PDB"] = String[]
    mmcif_dict["_atom_site.pdbx_PDB_model_num"] = String[]
    mmcif_dict["_atom_site.pdbx_PDB_ins_code"] = String[]

    # Dummy values to allow the convertion to BioStructures.MolecularStructure
    if molecular_structures
        mmcif_dict["_atom_site.label_alt_id"] = String[]
        mmcif_dict["_atom_site.pdbx_formal_charge"] = String[]
    end

    atom_id_counter = 1

    for res in residues
        for atom in res.atoms
            if molecular_structures || !label
                push!(mmcif_dict["_atom_site.auth_asym_id"], res.id.chain)
                push!(mmcif_dict["_atom_site.auth_comp_id"], res.id.name)
                push!(mmcif_dict["_atom_site.auth_atom_id"], atom.atom)
            end
            if label
                push!(mmcif_dict["_atom_site.label_asym_id"], res.id.chain)
                push!(mmcif_dict["_atom_site.label_comp_id"], res.id.name)
                push!(mmcif_dict["_atom_site.label_atom_id"], atom.atom)
            end
            push!(mmcif_dict["_atom_site.id"], string(atom_id_counter))
            push!(mmcif_dict["_atom_site.auth_seq_id"], _resnumber(res.id.number))
            push!(mmcif_dict["_atom_site.label_seq_id"], _resnumber(res.id.PDBe_number))
            push!(mmcif_dict["_atom_site.Cartn_x"], string(atom.coordinates.x))
            push!(mmcif_dict["_atom_site.Cartn_y"], string(atom.coordinates.y))
            push!(mmcif_dict["_atom_site.Cartn_z"], string(atom.coordinates.z))
            push!(mmcif_dict["_atom_site.occupancy"], string(atom.occupancy))
            push!(mmcif_dict["_atom_site.B_iso_or_equiv"], atom.B)
            push!(mmcif_dict["_atom_site.type_symbol"], atom.element)
            push!(mmcif_dict["_atom_site.group_PDB"], res.id.group)
            push!(mmcif_dict["_atom_site.pdbx_PDB_model_num"], res.id.model)
            push!(mmcif_dict["_atom_site.pdbx_PDB_ins_code"], _inscode(res))

            # Dummy values to allow the convertion to BioStructures.MolecularStructure
            if molecular_structures
                push!(mmcif_dict["_atom_site.label_alt_id"], ".")
                push!(mmcif_dict["_atom_site.pdbx_formal_charge"], "?")
            end
            atom_id_counter += 1
        end
    end

    return mmcif_dict
end

function Utils.print_file(
    io::IO,
    residues::AbstractVector{PDBResidue},
    format::Type{MMCIFFile};
    label::Bool = false,
)
    mmcif_dict = _pdbresidues_to_mmcifdict(residues, label = label)
    writemmcif(io, mmcif_dict)
end
