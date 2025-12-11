"""
`MMCIFFile <: FileFormat`

macromolecular Crystallographic Information File (mmCIF) format.
"""
struct MMCIFFile <: FileFormat end

_clean_string(s::String) = replace(s, "." => "", "?" => "")

function _parse_mmcif_to_pdbresidues(mmcif_dict::BioStructures.MMCIFDict, label::Bool)
    # Choose the correct prefix based on the label argument
    prefix = if label
        "_atom_site.label"
    else
        "_atom_site.auth"
    end
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
    alt_ids = mmcif_dict["_atom_site.label_alt_id"] # Alternative location ID
    formal_charges =
        haskey(mmcif_dict, "_atom_site.pdbx_formal_charge") ?
        mmcif_dict["_atom_site.pdbx_formal_charge"] : ["?" for _ = 1:length(atom_names)]
    # Formal charge: net integer charge assigned; default is "0" as the regex is [+-]?[0-9]+
    # https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_atom_site.pdbx_formal_charge.html
    # If the value is missing, it is represented as "?" in the mmCIF file. 
    # The `_clean_string` function will convert "?" to "", which we then print as "?" 
    # thanks to the `_pdbresidues_to_mmcifdict` function.
    # https://mmcif.wwpdb.org/docs/tutorials/mechanics/pdbx-mmcif-syntax.html
    ins_codes = mmcif_dict["_atom_site.pdbx_PDB_ins_code"] # Insertion codes

    residues = PDBResidue[]
    current_residue_id = ""
    current_residue = PDBResidue(PDBResidueIdentifier("", "", "", "", "", ""), PDBAtom[])

    for i = 1:length(auth_seq_ids)
        pdb_number = string(auth_seq_ids[i], _clean_string(ins_codes[i]))
        pdbe_number = _clean_string(label_seq_ids[i])
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
            _clean_string(alt_ids[i]),
            _clean_string(formal_charges[i]),
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
)::Vector{PDBResidue}
    mmcif_dict = BioStructures.MMCIFDict(io)

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

function _pdbresidues_to_mmcifdict(residues::Vector{PDBResidue}; label::Bool = false)
    n_atoms = sum(length(res.atoms) for res in residues)

    mmcif_dict = BioStructures.MMCIFDict()

    if label
        mmcif_dict["_atom_site.label_asym_id"] = Vector{String}(undef, n_atoms)
        mmcif_dict["_atom_site.label_comp_id"] = Vector{String}(undef, n_atoms)
        mmcif_dict["_atom_site.label_atom_id"] = Vector{String}(undef, n_atoms)
    else
        mmcif_dict["_atom_site.auth_asym_id"] = Vector{String}(undef, n_atoms)
        mmcif_dict["_atom_site.auth_comp_id"] = Vector{String}(undef, n_atoms)
        mmcif_dict["_atom_site.auth_atom_id"] = Vector{String}(undef, n_atoms)
    end
    mmcif_dict["_atom_site.id"] = Vector{String}(undef, n_atoms)
    mmcif_dict["_atom_site.auth_seq_id"] = Vector{String}(undef, n_atoms)
    mmcif_dict["_atom_site.label_seq_id"] = Vector{String}(undef, n_atoms)
    mmcif_dict["_atom_site.Cartn_x"] = Vector{String}(undef, n_atoms)
    mmcif_dict["_atom_site.Cartn_y"] = Vector{String}(undef, n_atoms)
    mmcif_dict["_atom_site.Cartn_z"] = Vector{String}(undef, n_atoms)
    mmcif_dict["_atom_site.occupancy"] = Vector{String}(undef, n_atoms)
    mmcif_dict["_atom_site.B_iso_or_equiv"] = Vector{String}(undef, n_atoms)
    mmcif_dict["_atom_site.type_symbol"] = Vector{String}(undef, n_atoms)
    mmcif_dict["_atom_site.group_PDB"] = Vector{String}(undef, n_atoms)
    mmcif_dict["_atom_site.pdbx_PDB_model_num"] = Vector{String}(undef, n_atoms)
    mmcif_dict["_atom_site.pdbx_PDB_ins_code"] = Vector{String}(undef, n_atoms)
    mmcif_dict["_atom_site.label_alt_id"] = Vector{String}(undef, n_atoms)
    mmcif_dict["_atom_site.pdbx_formal_charge"] = Vector{String}(undef, n_atoms)

    i = 1
    for res in residues
        for atom in res.atoms
            if label
                mmcif_dict["_atom_site.label_asym_id"][i] = res.id.chain
                mmcif_dict["_atom_site.label_comp_id"][i] = res.id.name
                mmcif_dict["_atom_site.label_atom_id"][i] = atom.atom
            else
                mmcif_dict["_atom_site.auth_asym_id"][i] = res.id.chain
                mmcif_dict["_atom_site.auth_comp_id"][i] = res.id.name
                mmcif_dict["_atom_site.auth_atom_id"][i] = atom.atom
            end
            mmcif_dict["_atom_site.id"][i] = string(i)
            mmcif_dict["_atom_site.auth_seq_id"][i] = _resnumber(res.id.number)
            mmcif_dict["_atom_site.label_seq_id"][i] = _resnumber(res.id.PDBe_number)
            mmcif_dict["_atom_site.Cartn_x"][i] = string(atom.coordinates.x)
            mmcif_dict["_atom_site.Cartn_y"][i] = string(atom.coordinates.y)
            mmcif_dict["_atom_site.Cartn_z"][i] = string(atom.coordinates.z)
            mmcif_dict["_atom_site.occupancy"][i] = string(atom.occupancy)
            mmcif_dict["_atom_site.B_iso_or_equiv"][i] = atom.B
            mmcif_dict["_atom_site.type_symbol"][i] = atom.element
            mmcif_dict["_atom_site.group_PDB"][i] = res.id.group
            mmcif_dict["_atom_site.pdbx_PDB_model_num"][i] = res.id.model
            mmcif_dict["_atom_site.pdbx_PDB_ins_code"][i] = _inscode(res)
            mmcif_dict["_atom_site.label_alt_id"][i] =
                isempty(atom.alt_id) ? "." : atom.alt_id
            mmcif_dict["_atom_site.pdbx_formal_charge"][i] =
                isempty(atom.charge) ? "?" : atom.charge
            i += 1
        end
    end

    return mmcif_dict
end

function Utils.print_file(
    io::IO,
    residues::AbstractVector{PDBResidue},
    format::Type{MMCIFFile};
    label::Bool = true,
)
    mmcif_dict = _pdbresidues_to_mmcifdict(residues, label = label)
    writemmcif(io, mmcif_dict)
end

"""
    BioStructures.MMCIFDict(residues::Vector{PDBResidue}; label::Bool = false)

Create a `BioStructures.MMCIFDict` from a vector of `PDBResidue`s. Set
`label = true` to fill the CIF dictionary using the `label_` fields instead of
the `auth_` fields.
"""
function BioStructures.MMCIFDict(residues::Vector{PDBResidue}; label::Bool = false)
    _pdbresidues_to_mmcifdict(residues; label = label)
end

"""
    Base.convert(::Type{Vector{PDBResidue}}, mmcif_dict::BioStructures.MMCIFDict)

Convert a `MMCIFDict` into a vector of `PDBResidue`s using the `label_` fields
present in the dictionary.
"""
function Base.convert(::Type{Vector{PDBResidue}}, mmcif_dict::BioStructures.MMCIFDict)
    label = haskey(mmcif_dict, "_atom_site.label_asym_id")
    _parse_mmcif_to_pdbresidues(mmcif_dict, label)
end

"""
    Base.convert(::Type{BioStructures.MMCIFDict}, residues::Vector{PDBResidue})

Return a `MMCIFDict` representation of `residues` using the `auth_` fields.
"""
Base.convert(::Type{BioStructures.MMCIFDict}, residues::Vector{PDBResidue}) =
    BioStructures.MMCIFDict(residues)
