# This file defines the functions to convert between MIToS.PDB and BioStructures types.

function BioStructures.MolecularStructure(residues::Vector{PDBResidue})
    mmcifdict = _pdbresidues_to_mmcifdict(residues, molecular_structures = true)
    MolecularStructure(mmcifdict)
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
        )
        push!(atoms, atom_obj)
    end

    PDBResidue(residue_id, atoms)
end
