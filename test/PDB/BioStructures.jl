using Test
using MIToS.PDB
using MIToS.Utils # For parse_file
using BioStructures
import MIToS.PDB: _pdbresidues_to_mmcifdict # Internal function for reference

# Define this to make BioStructures constructor accessible for testing
# This might be needed if the constructor isn't automatically resolved
# to the one in MIToS.PDB's scope when called with PDBResidues.
# However, Julia's method extension should handle this:
# The `BioStructures.MolecularStructure(residues::Vector{PDBResidue})` in
# `src/PDB/BioStructures.jl` extends the BioStructures.MolecularStructure function.

# Helper function to compare two MolecularStructure objects
function are_molecular_structures_equal(s1::MolecularStructure, s2::MolecularStructure; atol::Float64 = 1e-5)
    # Compare number of models
    if length(models(s1)) != length(models(s2))
        @warn "Different number of models: $(length(models(s1))) vs $(length(models(s2)))"
        return false
    end

    m1_dict = models(s1)
    m2_dict = models(s2)

    # Iterate through models - assuming model IDs should match and be in order
    # BioStructures models are typically OrderedDict, so order should be preserved
    model_ids1 = collect(keys(m1_dict))
    model_ids2 = collect(keys(m2_dict))

    if model_ids1 != model_ids2
        @warn "Model IDs differ or are in different order: $model_ids1 vs $model_ids2"
        return false
    end

    for model_id in model_ids1
        model1 = m1_dict[model_id]
        model2 = m2_dict[model_id]

        if modelnumber(model1) != modelnumber(model2) # Should be same as model_id
            @warn "Model numbers differ for model ID $model_id: $(modelnumber(model1)) vs $(modelnumber(model2))"
            return false
        end

        # Compare number of chains
        if length(chains(model1)) != length(chains(model2))
            @warn "Different number of chains in model $model_id: $(length(chains(model1))) vs $(length(chains(model2)))"
            return false
        end

        c1_dict = chains(model1)
        c2_dict = chains(model2)
        chain_ids1 = collect(keys(c1_dict))
        chain_ids2 = collect(keys(c2_dict))

        if chain_ids1 != chain_ids2
            @warn "Chain IDs differ or are in different order in model $model_id: $chain_ids1 vs $chain_ids2"
            return false
        end

        for chain_id in chain_ids1
            chain1 = c1_dict[chain_id]
            chain2 = c2_dict[chain_id]

            if chainid(chain1) != chainid(chain2) # Should be same as chain_id
                @warn "Chain IDs differ for chain $chain_id in model $model_id: $(chainid(chain1)) vs $(chainid(chain2))"
                return false
            end

            # Compare number of residues
            # residues(chain) returns OrderedDict{ResidueIdentifier, AbstractResidue}
            res1_dict = BioStructures.residues(chain1) # Use BioStructures.residues to be explicit
            res2_dict = BioStructures.residues(chain2)
            
            if length(res1_dict) != length(res2_dict)
                @warn "Different number of residues in chain $chain_id, model $model_id: $(length(res1_dict)) vs $(length(res2_dict))"
                return false
            end

            res_ids1 = collect(keys(res1_dict))
            res_ids2 = collect(keys(res2_dict))

            if res_ids1 != res_ids2 # ResidueIdentifier should be comparable
                @warn "Residue Identifiers differ or are in different order in chain $chain_id, model $model_id: $res_ids1 vs $res_ids2"
                return false
            end
            
            for res_id in res_ids1
                res1 = res1_dict[res_id]
                res2 = res2_dict[res_id]

                if resname(res1) != resname(res2)
                    @warn "Residue names differ for $res_id: $(resname(res1)) vs $(resname(res2))"
                    return false
                end
                if resnumber(res1) != resnumber(res2) # From ResidueIdentifier
                    @warn "Residue numbers differ for $res_id: $(resnumber(res1)) vs $(resnumber(res2))"
                    return false
                end
                if inscode(res1) != inscode(res2) # From ResidueIdentifier
                    @warn "Insertion codes differ for $res_id: $(inscode(res1)) vs $(inscode(res2))"
                    return false
                end
                if ishetero(res1) != ishetero(res2)
                    @warn "Hetero flags differ for $res_id: $(ishetero(res1)) vs $(ishetero(res2))"
                    return false
                end

                # Compare atoms - collectatoms handles disordered residues
                atoms1 = collectatoms(res1) # Returns Vector{AbstractAtom}
                atoms2 = collectatoms(res2)

                if length(atoms1) != length(atoms2)
                    @warn "Different number of atoms in residue $res_id, chain $chain_id, model $model_id: $(length(atoms1)) vs $(length(atoms2))"
                    return false
                end

                for (i, (atom1, atom2)) in enumerate(zip(atoms1, atoms2))
                    if atomname(atom1) != atomname(atom2)
                        @warn "Atom names differ at index $i for $res_id: $(atomname(atom1)) vs $(atomname(atom2))"
                        return false
                    end
                    if !isapprox(BioStructures.x(atom1), BioStructures.x(atom2); atol=atol) ||
                       !isapprox(BioStructures.y(atom1), BioStructures.y(atom2); atol=atol) ||
                       !isapprox(BioStructures.z(atom1), BioStructures.z(atom2); atol=atol)
                        @warn "Atom coordinates differ for $(atomname(atom1)) in $res_id: ($(BioStructures.x(atom1)),$(BioStructures.y(atom1)),$(BioStructures.z(atom1))) vs ($(BioStructures.x(atom2)),$(BioStructures.y(atom2)),$(BioStructures.z(atom2)))"
                        return false
                    end
                    if !isapprox(occupancy(atom1), occupancy(atom2); atol=atol)
                        @warn "Atom occupancies differ for $(atomname(atom1)) in $res_id: $(occupancy(atom1)) vs $(occupancy(atom2))"
                        return false
                    end
                    if !isapprox(tempfactor(atom1), tempfactor(atom2); atol=atol)
                        @warn "Atom tempfactors differ for $(atomname(atom1)) in $res_id: $(tempfactor(atom1)) vs $(tempfactor(atom2))"
                        return false
                    end
                    if altlocid(atom1) != altlocid(atom2)
                        @warn "Atom altlocids differ for $(atomname(atom1)) in $res_id: $(altlocid(atom1)) vs $(altlocid(atom2))"
                        return false
                    end
                    if charge(atom1) != charge(atom2)
                        @warn "Atom charges differ for $(atomname(atom1)) in $res_id: $(charge(atom1)) vs $(charge(atom2))"
                        return false
                    end
                    if BioStructures.element(atom1) != BioStructures.element(atom2)
                        @warn "Atom elements differ for $(atomname(atom1)) in $res_id: $(BioStructures.element(atom1)) vs $(BioStructures.element(atom2))"
                        return false
                    end
                    # Note: is_hetero_atom field in BioStructures.Atom is not directly checked here
                    # as it's often inferred or used internally. The residue's ishetero() check is primary.
                end
            end
        end
    end
    return true # All checks passed
end


@testset "Direct PDBResidue to MolecularStructure Conversion" begin
    # Path to the test PDB file, relative to the test/PDB directory
    # Assuming this test file is in test/PDB/ and data is in test/data/
    pdb_filepath = "../data/1AS5.pdb" 

    if !isfile(pdb_filepath)
        @warn "Test PDB file not found at $pdb_filepath. Skipping tests."
        # To ensure tests don't error out if file is missing, but still report.
        # A better way might be to download it if not present, or have it as a fixture.
        @test_broken false # Indicates that the test setup failed.
        return # Exit testset
    end

    mitos_pdb_residues = MIToS.Utils.parse_file(pdb_filepath, PDBFile)

    # 1. Create MolecularStructure using the new direct constructor
    # The constructor being tested is MIToS.PDB.BioStructures.MolecularStructure,
    # which extends ::BioStructures.MolecularStructure.
    # We can pass a structure_id for consistency if needed, e.g., "TEST_ID"
    # The default in the implementation is "MIToSStructure"
    struct_new = BioStructures.MolecularStructure(mitos_pdb_residues; structure_id = "1AS5_NEW")

    # 2. Create a reference MolecularStructure using the old mmCIF-based method
    # This uses the BioStructures.MolecularStructure(mmcif_dict) constructor from BioStructures.jl
    mmcif_dict = _pdbresidues_to_mmcifdict(mitos_pdb_residues, molecular_structures = true)
    
    # The BioStructures.MolecularStructure constructor from mmcif_dict might derive its ID
    # from the dict. For this test, we are not comparing structure IDs directly in
    # are_molecular_structures_equal, so this is fine.
    # If mmcif_dict contains an ID, BioStructures will use it.
    # If not, it might default to something like "default".
    struct_old_ref = BioStructures.MolecularStructure(mmcif_dict) 
    # To make it more comparable if its ID matters, one might need to inspect mmcif_dict
    # and pass that ID, or ensure BioStructures.MolecularStructure(mmcif_dict, "ID") exists.
    # For now, are_molecular_structures_equal will focus on content.

    @test are_molecular_structures_equal(struct_new, struct_old_ref) "Molecular structures created by new and old methods are not equal."
end
