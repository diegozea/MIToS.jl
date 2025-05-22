using MIToS.PDB
using MIToS.Utils # For All if needed
using RecipesBase
using Test
using StaticArrays # For Coordinates if creating them directly

@testset "PDB Plot Recipes" begin
    # Helper to create a dummy PDBResidue vector
    function create_dummy_pdb_residues(nres)
        residues = PDBResidue[]
        for i in 1:nres
            # Create a C-alpha atom. Coordinates are Float64 as per PDBAtom constructor.
            atom_ca = PDBAtom(Coordinates(Float64(i), 2.0*i, 3.0*i), "CA", "C", 1.0, "0.0", "", "")
            # Create a Nitrogen atom
            atom_n = PDBAtom(Coordinates(0.5*i, 1.5*i, 2.5*i), "N", "N", 1.0, "0.0", "", "")
            
            # Alternate chain and group
            chain_id = string('A' + mod(i-1, 2))
            # Mix ATOM and HETATM: Res 1 (ATOM), Res 2 (HETATM), Res 3 (ATOM), Res 4 (HETATM), Res 5 (ATOM)
            group_type = mod(i-1, 2) == 0 ? "ATOM" : "HETATM"
            
            residue_id = PDBResidueIdentifier("", string(i), "ALA", group_type, "", chain_id) # Inscode is empty string
            
            # For some residues, make CA atom missing to test NaN handling
            # Residue 3 will have no CA atom
            atoms_list = if i == 3 
                PDBAtom[atom_n] # Only N atom, no CA
            else
                PDBAtom[atom_ca, atom_n] # Both CA and N atoms
            end
            push!(residues, PDBResidue(residue_id, atoms_list))
        end
        return residues
    end

    @testset "PDBResidue Vector Plot Recipe" begin
        residues = create_dummy_pdb_residues(5) # Create 5 residues
        # Expected: Res1 (ATOM,A), Res2 (HETATM,B), Res3 (ATOM,A, no CA), Res4 (HETATM,B), Res5 (ATOM,A)
        
        plot_data_list = RecipesBase.apply_recipe(Dict{Symbol,Any}(), residues)
        @test !isempty(plot_data_list)
        @test plot_data_list isa Vector{RecipesBase.PlotData}
        plot_data = plot_data_list[1]
        attrs = plot_data.plotattributes

        # Test :group attribute: Should only include chains from "ATOM" residues
        # Expected ATOM residues: 1 (Chain A), 3 (Chain A), 5 (Chain A)
        expected_chains = [r.id.chain for r in residues if r.id.group == "ATOM"]
        @test attrs[:group] == expected_chains
        @test attrs[:group] == ["A", "A", "A"] # Explicitly for this setup

        # Test data arguments (x, y, z coordinates from CAmatrix)
        ca_matrix = CAmatrix(residues) # Expected to have 5 rows, 3 columns
                                       # Row 3 should be NaNs due to missing CA in residue 3
        
        @test size(ca_matrix) == (5, 3) # Verify shape of CAmatrix output
        
        # Check coordinates
        @test plot_data.args[1] == ca_matrix[:, 1] # X coordinates
        @test plot_data.args[2] == ca_matrix[:, 2] # Y coordinates
        @test plot_data.args[3] == ca_matrix[:, 3] # Z coordinates
        
        # Check for NaNs if CA was missing (for residue 3 in helper)
        @test isnan(plot_data.args[1][3]) # X for residue 3
        @test isnan(plot_data.args[2][3]) # Y for residue 3
        @test isnan(plot_data.args[3][3]) # Z for residue 3

        # Check non-NaN values for other residues to be sure
        @test plot_data.args[1][1] == 1.0
        @test plot_data.args[2][1] == 2.0
        @test plot_data.args[3][1] == 3.0
        
        @test plot_data.args[1][2] == 2.0 # HETATM, should still have coords in CAmatrix
        @test plot_data.args[2][2] == 4.0
        @test plot_data.args[3][2] == 6.0
    end

    @testset "PDBResidue Vector Plot Recipe - Empty Input" begin
        residues_empty = PDBResidue[]
        plot_data_list_empty = RecipesBase.apply_recipe(Dict{Symbol,Any}(), residues_empty)
        @test !isempty(plot_data_list_empty)
        plot_data_empty = plot_data_list_empty[1]
        
        @test isempty(plot_data_empty.plotattributes[:group])
        @test isempty(plot_data_empty.args[1]) # X
        @test isempty(plot_data_empty.args[2]) # Y
        @test isempty(plot_data_empty.args[3]) # Z
    end
    
    @testset "PDBResidue Vector Plot Recipe - All HETATM" begin
        residues_het = PDBResidue[]
        for i in 1:3
            atom_ca = PDBAtom(Coordinates(Float64(i), 2.0*i, 3.0*i), "CA", "C", 1.0, "0.0", "", "")
            residue_id = PDBResidueIdentifier("", string(i), "ALA", "HETATM", "", "A") # Inscode empty
            push!(residues_het, PDBResidue(residue_id, PDBAtom[atom_ca]))
        end
        plot_data_list_het = RecipesBase.apply_recipe(Dict{Symbol,Any}(), residues_het)
        @test !isempty(plot_data_list_het)
        plot_data_het = plot_data_list_het[1]
        
        @test isempty(plot_data_het.plotattributes[:group]) # No "ATOM" groups, so group should be empty
        
        # Data should still be there as CAmatrix processes HETATM too if they have CA
        ca_matrix_het = CAmatrix(residues_het)
        @test plot_data_het.args[1] == ca_matrix_het[:, 1]
        @test plot_data_het.args[2] == ca_matrix_het[:, 2]
        @test plot_data_het.args[3] == ca_matrix_het[:, 3]
        @test length(plot_data_het.args[1]) == 3
    end
end
