# This file defines the functions to convert between MIToS.PDB and BioStructures types.

# Helper function to parse MIToS residue number string (e.g., "10A")
# into residue number (Int) and insertion code (Char).
function _get_resnum_inscode_from_mitos(mitos_res_number::String)
    if isempty(mitos_res_number)
        return (0, ' ') # Should not happen with valid PDBResidue
    end
    last_char = mitos_res_number[end]
    if isletter(last_char)
        resnum_str = mitos_res_number[1:end-1]
        inscode = last_char
    else
        resnum_str = mitos_res_number
        inscode = ' '
    end
    return (parse(Int, resnum_str), inscode)
end

# Helper function to parse charge string to Int, defaulting to 0.
function _parse_charge(charge_str::String)
    s = strip(charge_str)
    if isempty(s)
        return 0
    end
    # Handle forms like "1+", "+1", "1", "-1", "1-"
    val_str = filter(isdigit, s)
    if isempty(val_str)
        return 0 # No digits found
    end
    val = parse(Int, val_str)
    if occursin('-', s)
        return -val
    else
        return val
    end
end

function BioStructures.MolecularStructure(residues::Vector{PDBResidue}; structure_id::String="MIToSStructure")
    models_dict = OrderedDict{String,Model}()

    for mitos_res in residues
        model_id_str = mitos_res.id.model
        chain_id_str = mitos_res.id.chain
        
        resnum, inscode = _get_resnum_inscode_from_mitos(mitos_res.id.number)
        res_name_str = mitos_res.id.name
        is_hetero_residue = mitos_res.id.group == "HETATM"

        # Get or create Model
        model = get!(models_dict, model_id_str) do
            Model(model_id_str, OrderedDict{String,Chain}())
        end

        # Get or create Chain
        chain = get!(model.chains, chain_id_str) do
            Chain(chain_id_str, OrderedDict{ResidueIdentifier,AbstractResidue}())
        end

        # Create Atom list for the current residue
        atom_list = Vector{BioStructures.Atom}()
        for mitos_atom in mitos_res.atoms
            coords = Coordinates(mitos_atom.coor.x, mitos_atom.coor.y, mitos_atom.coor.z)
            atom_name_str = mitos_atom.atom
            # alt_loc_id must be Char, default to ' ' if empty
            alt_loc_id_char = isempty(mitos_atom.alt_id) ? ' ' : first(mitos_atom.alt_id)
            charge_int = _parse_charge(mitos_atom.charge)
            
            b_factor_float = 0.0 # Default B-factor
            try
                b_factor_float = parse(Float64, mitos_atom.B)
            catch e
                # Add warning or error handling if necessary for malformed B-factor
                # For now, defaults to 0.0 if parsing fails
            end

            occupancy_float = mitos_atom.occupancy # Already Float64 in MIToS.PDBAtom
            element_str = mitos_atom.element
            
            # Atom's hetero status is typically inherited from the residue
            is_hetero_atom = is_hetero_residue 

            bs_atom = BioStructures.Atom(
                coords,
                atom_name_str,
                alt_loc_id_char,
                charge_int,
                b_factor_float,
                occupancy_float,
                element_str,
                is_hetero_atom # BioStructures.Atom uses this field
            )
            push!(atom_list, bs_atom)
        end

        # Create ResidueIdentifier for BioStructures
        # Based on BioStructures.jl, ResidueIdentifier is typically (resnumber, inscode)
        # resname is stored in the Residue object itself.
        bs_res_id = ResidueIdentifier(resnum, inscode)

        # Create BioStructures.Residue
        # TODO: Handle DisorderedResidue if necessary. For now, creating Residue.
        bs_residue = Residue(
            bs_res_id,
            res_name_str,
            atom_list,
            is_hetero_residue
        )
        
        # Add residue to chain
        chain.residues[bs_res_id] = bs_residue
    end

    return MolecularStructure(structure_id, models_dict)
end

function Base.convert(::Type{Vector{PDBResidue}}, struc::MolecularStructure)
    _molecularstructure_to_pdbresidues(struc)
end

function _molecularstructure_to_pdbresidues(struc::MolecularStructure)
    vector_res = PDBResidue[]
    for model in values(models(struc))
        for chain in values(chains(model))
            for res in values(BioStructures.residues(chain))
                # Check if the residue is a DisorderedResidue
                if isdisorderedres(res)
                    for res_name in resnames(res)
                        dis_res = disorderedres(res, res_name)
                        push!(vector_res, _create_pdbresidue(dis_res, model, chain))
                    end
                else
                    push!(vector_res, _create_pdbresidue(res, model, chain))
                end
            end
        end
    end
    vector_res
end

function _create_pdbresidue(res, model::Model, chain::Chain)
    residue_id = PDBResidueIdentifier(
        "",  # PDBe_number not available
        string(resnumber(res)) * inscode(res),
        resname(res),
        ishetero(res) ? "HETATM" : "ATOM",
        string(modelnumber(model)),
        chainid(chain),
    )

    atoms = PDBAtom[]

    for atom in collectatoms(res)
        atom_obj = PDBAtom(
            Coordinates(
                BioStructures.x(atom),
                BioStructures.y(atom),
                BioStructures.z(atom),
            ),
            atomname(atom),
            element(atom),
            occupancy(atom),
            string(tempfactor(atom)),
            string(altlocid(atom)),
            strip(charge(atom, strip = false)),
        )
        push!(atoms, atom_obj)
    end

    PDBResidue(residue_id, atoms)
end
