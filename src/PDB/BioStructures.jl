# This file defines the functions to convert between MIToS.PDB and BioStructures types.

"""
    BioStructures.MolecularStructure(residues::Vector{PDBResidue})

Build a `BioStructures.MolecularStructure` from a vector of `PDBResidue`
objects without creating an intermediate mmCIF dictionary.
"""
function BioStructures.MolecularStructure(residues::Vector{PDBResidue})
    struc = MolecularStructure()
    serial = 1
    for res in residues
        model_n = parse(Int, res.id.model)
        if !haskey(models(struc), model_n)
            struc[model_n] = Model(model_n, struc)
        end
        mo = struc[model_n]
        chain_id = res.id.chain
        if !haskey(chains(mo), chain_id)
            mo[chain_id] = Chain(chain_id, mo)
        end
        res_number = parse(Int, replace(res.id.number, r"[A-Za-z]" => ""))
        ins_match = match(r"[A-Za-z]$", res.id.number)
        ins_code = ins_match === nothing ? ' ' : ins_match.match[1]
        hetero = res.id.group != "ATOM"
        for at in res.atoms
            coords = [at.coordinates.x, at.coordinates.y, at.coordinates.z]
            alt_id = isempty(at.alt_id) ? '\0' : at.alt_id[1]
            b = tryparse(Float64, at.B)
            temp = b === nothing ? 0.0 : b
            rec = AtomRecord(
                hetero,
                serial,
                at.atom,
                alt_id,
                res.id.name,
                chain_id,
                res_number,
                ins_code,
                coords,
                at.occupancy,
                temp,
                at.element,
                at.charge,
            )
            BioStructures.unsafe_addatomtomodel!(mo, rec, remove_disorder = false)
            serial += 1
        end
    end
    BioStructures.fixlists!(struc)
    return struc
end

"""
    Base.convert(::Type{Vector{PDBResidue}}, struc::MolecularStructure)

Return a vector of `PDBResidue`s representing the given
`MolecularStructure`. All models, chains, residues and atoms are
retained and alternate residues or atoms are expanded.
"""
function Base.convert(::Type{Vector{PDBResidue}}, struc::MolecularStructure)
    _molecularstructure_to_pdbresidues(struc)
end

function _molecularstructure_to_pdbresidues(struc::MolecularStructure)
    nres = BioStructures.countresidues(struc; expand_disordered = true)
    residues = Vector{PDBResidue}()
    sizehint!(residues, nres)

    for model in values(BioStructures.models(struc))
        model_id = string(BioStructures.modelnumber(model))
        for chain in values(BioStructures.chains(model))
            chain_id = BioStructures.chainid(chain)
            for logical_res in values(BioStructures.residues(chain))
                if BioStructures.isdisorderedres(logical_res)
                    for rname in BioStructures.resnames(logical_res)
                        res = BioStructures.disorderedres(logical_res, rname)
                        push!(residues, _create_pdbresidue(res, model_id, chain_id))
                    end
                else
                    push!(residues, _create_pdbresidue(logical_res, model_id, chain_id))
                end
            end
        end
    end

    return residues
end

@inline function _create_pdbresidue(res, model_id::String, chain_id::String)
    residue_id = PDBResidueIdentifier(
        "",
        string(BioStructures.resnumber(res)) *
        (BioStructures.inscode(res) == ' ' ? "" : string(BioStructures.inscode(res))),
        BioStructures.resname(res),
        BioStructures.ishetero(res) ? "HETATM" : "ATOM",
        model_id,
        chain_id,
    )

    atoms = PDBAtom[]
    @inbounds for at in res
        _push_atom!(atoms, at)
        if BioStructures.isdisorderedatom(at)
            for alt_id in BioStructures.altlocids(at)
                alt_id == BioStructures.altlocid(at) && continue
                altat = res[BioStructures.atomname(at)][alt_id]
                _push_atom!(atoms, altat)
            end
        end
    end

    return PDBResidue(residue_id, atoms)
end

@inline function _push_atom!(atoms::Vector{PDBAtom}, at)
    push!(
        atoms,
        PDBAtom(
            Coordinates(BioStructures.x(at), BioStructures.y(at), BioStructures.z(at)),
            BioStructures.atomname(at),
            BioStructures.element(at),
            BioStructures.occupancy(at),
            string(BioStructures.tempfactor(at)),
            BioStructures.altlocid(at) == '\0' ? "" : string(BioStructures.altlocid(at)),
            BioStructures.charge(at, strip = false),
        ),
    )
end
