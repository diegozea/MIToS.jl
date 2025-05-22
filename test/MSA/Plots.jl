using MIToS.MSA
using MIToS.Utils # For All if needed, or other utilities
using RecipesBase
using Test
using OrderedCollections # For NamedArray construction
using MIToS.Residues # Required for res"..."

@testset "MSA Plot Recipes" begin
    # Helper to create a dummy MSA for testing
    function create_dummy_msa(nseq, ncol)
        # Ensure we use valid residues for the MSA
        valid_residues = res"ARNDCQEGHILKMFPSTWYV-"
        data = NamedArray(rand(valid_residues, nseq, ncol)) # Use NamedArray directly
        seqnames = ["Seq" * string(i) for i in 1:nseq]
        colnames = [string(j) for j in 1:ncol] # Corrected to use j for columns

        # Set names for the NamedArray
        setnames!(data, seqnames, 1)
        setnames!(data, colnames, 2)
        
        # For basic recipe testing, MultipleSequenceAlignment should suffice
        return MultipleSequenceAlignment(data)
    end

    # Expected _residue_labels for verification
    # Based on src/MSA/Plots.jl: map(string, reverse!(res"ARNDCQEGHILKMFPSTWYV-"))
    # It's an internal detail, so replicating its expected value is more robust than accessing it directly if not exported.
    expected_residue_labels = map(string, Residue['-'], Residue['V'], Residue['Y'], Residue['W'], Residue['T'], Residue['S'], Residue['P'], Residue['F'], Residue['M'], Residue['K'], Residue['L'], Residue['I'], Residue['H'], Residue['G'], Residue['E'], Residue['Q'], Residue['C'], Residue['D'], Residue['N'], Residue['R'], Residue['A'])


    @testset "Basic Heatmap Recipe" begin
        msa = create_dummy_msa(5, 10) # Small MSA
        plot_data_list = RecipesBase.apply_recipe(Dict{Symbol,Any}(), msa)
        @test !isempty(plot_data_list)
        @test plot_data_list isa Vector{RecipesBase.PlotData}
        plot_data = plot_data_list[1]
        attrs = plot_data.plotattributes

        @test attrs[:seriestype] == :heatmap
        @test attrs[:yflip] == true
        @test attrs[:grid] == false
        @test attrs[:foreground_color_border] == nothing
        @test attrs[:foreground_color_axis] == nothing
        @test attrs[:linewidth] == 0
        @test issetequal(attrs[:zdiscrete_values], expected_residue_labels) # Use issetequal for order-agnostic comparison

        # Test data arguments
        @test plot_data.args[1] == (1:ncolumns(msa)) # Corrected: ensure it's a tuple for range
        @test plot_data.args[2] == sequencenames(msa)
        # The third argument to the heatmap should be the matrix of residues, 
        # but the recipe converts them to strings first.
        # The recipe itself passes `getresidues(msa)` and `RecipesPipeline.jl` handles conversion based on zdiscrete_values
        # However, the provided example tests for `map(string, getresidues(msa))`
        # Let's check what the recipe actually produces for args[3] after applying the recipe.
        # The recipe passes `getresidues(msa)` to the heatmap series.
        # The `zdiscrete_values` handles the mapping from Residue to color.
        # So, args[3] should be the actual residue matrix.
        @test plot_data.args[3] == getresidues(msa)
        
        # Test yticks for small MSA (should not be set by recipe)
        @test !haskey(attrs, :yticks)
        @test !haskey(attrs, :html_output_format) # html_output_format should not be set for small MSAs
    end

    @testset "Y-axis Ticks for Large MSA" begin
        msa_large = create_dummy_msa(25, 10) # Large MSA (e.g., > 20 sequences)
        plot_data_list_large = RecipesBase.apply_recipe(Dict{Symbol,Any}(), msa_large)
        @test !isempty(plot_data_list_large)
        plot_data_large = plot_data_list_large[1]
        attrs_large = plot_data_large.plotattributes

        @test haskey(attrs_large, :yticks)
        nseq_large = nsequences(msa_large)
        # As per Plots.jl recipe logic: step = div(nseq_large, 20)
        step_large = div(nseq_large, 20) 
        expected_yticks_ticks = 1:step_large:nseq_large
        expected_yticks_labels = sequencenames(msa_large)[expected_yticks_ticks] # Ticks should align with labels

        @test attrs_large[:yticks][1] == expected_yticks_ticks
        @test attrs_large[:yticks][2] == expected_yticks_labels
        
        @test attrs_large[:html_output_format] == :png
    end

    @testset "AnnotatedMultipleSequenceAlignment Recipe" begin
        # Test with AnnotatedMultipleSequenceAlignment to ensure it also works
        annot_msa_data = rand(Residue, 3, 5)
        annot_seqnames = ["AnnSeq" * string(i) for i in 1:3]
        annot_colnames = [string(j) for j in 1:5]
        
        msa_part = MultipleSequenceAlignment(NamedArray(annot_msa_data, 
                                                        (OrderedDict(k => i for (i,k) in enumerate(annot_seqnames)), 
                                                         OrderedDict(k => i for (i,k) in enumerate(annot_colnames)))))
        
        # Create dummy annotations
        file_annot = OrderedDict("Version" => "1.0")
        col_annot = MSAAnnotations.Annotations("ColAnnot", OrderedDict(string(i) => string("Col", i) for i in 1:5))
        seq_annot = MSAAnnotations.Annotations("SeqAnnot", OrderedDict("AnnSeq1" => "Annotation S1"))
        
        annot_msa = AnnotatedMultipleSequenceAlignment(msa_part,
                                                       file_annot,
                                                       OrderedDict("Col" => col_annot), # Per Column Annotations
                                                       OrderedDict("Seq" => seq_annot)) # Per Sequence Annotations

        plot_data_list_annot = RecipesBase.apply_recipe(Dict{Symbol,Any}(), annot_msa)
        @test !isempty(plot_data_list_annot)
        plot_data_annot = plot_data_list_annot[1]
        attrs_annot = plot_data_annot.plotattributes

        @test attrs_annot[:seriestype] == :heatmap
        @test attrs_annot[:yflip] == true
        # Check a few key attributes to ensure the recipe was applied
        @test plot_data_annot.args[1] == (1:ncolumns(annot_msa))
        @test plot_data_annot.args[2] == sequencenames(annot_msa)
        @test plot_data_annot.args[3] == getresidues(annot_msa)
    end
end
