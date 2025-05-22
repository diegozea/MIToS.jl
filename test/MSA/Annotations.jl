using MIToS.MSA # Imports Annotations, MSAAnnotations
using Test
using OrderedCollections # For OrderedDict

@testset "Annotations Show/Print" begin

    @testset "Empty Annotations" begin
        ann = MSAAnnotations()
        # sprint(show, ...) for a type that defines show(io, mime, x) and show(io, x)
        # will use the 2-arg version if no MIME is specified.
        # For MSAAnnotations, show(io, ann) calls print_annotations(io, ann).
        
        str_output_show = sprint(show, ann)
        # Based on print_annotations, if all fields are empty, nothing is printed.
        @test str_output_show == ""

        str_output_print = sprint(print, ann)
        @test str_output_print == ""
    end

    @testset "Annotations with Content" begin
        # Define annotations similar to Stockholm format examples
        file_annot = OrderedDict{String,Vector{String}}(
            "ID" => ["TestPfam_ID"],
            "AC" => ["PF00000.1"]
        )
        # Sequence annotations: Key is sequence ID, value is Dict of feature -> values
        seq_annot = OrderedDict{String,OrderedDict{String,Vector{String}}}(
            "SEQ1/1-10" => OrderedDict("DR" => ["PDB;", "1XYZ A;"], "OS" => ["Human"]),
            "SEQ2/2-12"  => OrderedDict("DR" => ["GO;", "GO:1234567;"])
        )
        # Column annotations: Key is column ID, value is Dict of feature -> values
        # The prompt example "#=GC SS_cons HHH..EEE" implies a different structure for GC printing.
        # However, MSAAnnotations stores GC per column ID.
        # print_annotations iterates through `ann.col` which is `OrderedDict{String, OrderedDict{String, Vector{String}}}`
        # So, GC lines will be like #=GC <col_id> <feature> <value...>
        # To get something like "#=GC SS_cons HHHHEEE", the MSAAnnotations.col would need to be:
        # OrderedDict("SS_cons" => ["H","H","H","H","E","E","E"])
        # and print_annotations would need special handling for such keys.
        # The current print_annotations doesn't do this; it iterates `ann.col` as defined.
        # Let's test the actual behavior of print_annotations.
        col_annot = OrderedDict{String,OrderedDict{String,Vector{String}}}(
            "Col_1" => OrderedDict("FeatureA" => ["Val1"], "FeatureB" => ["ValX"]),
            "Col_3" => OrderedDict("FeatureA" => ["Val2"])
        )
        # Residue annotations: Key is (seq_id, col_id) tuple, value is Dict of feature -> values
        res_annot = OrderedDict{Tuple{String,String},OrderedDict{String,Vector{String}}}(
            ("SEQ1/1-10", "Col_1") => OrderedDict("Interaction" => ["SEQ2/2-12:Col_2"], "Charge" => ["+1"]),
            ("SEQ2/2-12", "Col_3") => OrderedDict("Modification" => ["Phosphorylation"])
        )

        ann = MSAAnnotations(file_annot, seq_annot, col_annot, res_annot)
        
        # sprint(show, ann) will call print_annotations(stdout, ann) internally.
        # We are testing the output of print_annotations.
        str_output = sprint(show, ann)

        # Expected File annotations
        @test occursin("#=GF ID TestPfam_ID", str_output)
        @test occursin("#=GF AC PF00000.1", str_output)
        
        # Expected Sequence annotations
        @test occursin("#=GS SEQ1/1-10 DR PDB; 1XYZ A;", str_output) # Values joined by space
        @test occursin("#=GS SEQ1/1-10 OS Human", str_output)
        @test occursin("#=GS SEQ2/2-12 DR GO; GO:1234567;", str_output)
        
        # Expected Column annotations (as per current print_annotations logic)
        @test occursin("#=GC Col_1 FeatureA Val1", str_output)
        @test occursin("#=GC Col_1 FeatureB ValX", str_output)
        @test occursin("#=GC Col_3 FeatureA Val2", str_output)
        
        # Expected Residue annotations
        @test occursin("#=GR SEQ1/1-10 Col_1 Interaction SEQ2/2-12:Col_2", str_output)
        @test occursin("#=GR SEQ1/1-10 Col_1 Charge +1", str_output)
        @test occursin("#=GR SEQ2/2-12 Col_3 Modification Phosphorylation", str_output)

        # Test that multiple values for a single annotation key are on the same line, space-separated
        file_annot_multi_val = OrderedDict("KW" => ["Keyword A", "Keyword B", "Keyword C"])
        ann_multi_gf = MSAAnnotations(file_annot_multi_val) # Only GF annotations
        str_output_multi_gf = sprint(show, ann_multi_gf)
        @test occursin("#=GF KW Keyword A Keyword B Keyword C", str_output_multi_gf)

        seq_annot_multi_val = OrderedDict("SEQ_MULTI" => OrderedDict("Features" => ["Feat1", "Feat2"]))
        ann_multi_gs = MSAAnnotations(OrderedDict{String,Vector{String}}(), seq_annot_multi_val) # Only GS
        str_output_multi_gs = sprint(show, ann_multi_gs)
        @test occursin("#=GS SEQ_MULTI Features Feat1 Feat2", str_output_multi_gs)
    end
    
    @testset "Annotations with only one type of content" begin
        # Test only GR annotations
        res_only_annot = OrderedDict{Tuple{String,String},OrderedDict{String,Vector{String}}}(
            ("S1", "C1") => OrderedDict("ResFeat" => ["Active"])
        )
        ann_gr_only = MSAAnnotations(
            OrderedDict{String,Vector{String}}(), # Empty GF
            OrderedDict{String,OrderedDict{String,Vector{String}}}(), # Empty GS
            OrderedDict{String,OrderedDict{String,Vector{String}}}(), # Empty GC
            res_only_annot # Only GR
        )
        str_output_gr_only = sprint(show, ann_gr_only)
        @test occursin("#=GR S1 C1 ResFeat Active", str_output_gr_only)
        @test !occursin("#=GF", str_output_gr_only)
        @test !occursin("#=GS", str_output_gr_only)
        @test !occursin("#=GC", str_output_gr_only)
    end
end
