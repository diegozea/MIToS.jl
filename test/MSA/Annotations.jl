@testset "Annotations" begin

    @testset "Empty annotations" begin

        annot = Annotations()

        @test length(annot.file) == 0
        @test length(annot.sequences) == 0
        @test length(annot.columns) == 0
        @test length(annot.residues) == 0

        @test length(annot) == 0

        @test ncolumns(annot) == -1

        @test isempty(annot)
    end

    @testset "Getters & Setters" begin

        annot = Annotations()
        example_str = "CCCHHHHHHHHHHHHHHHEEEEEEEEEEEEEEEEHHH"

        setannotresidue!(annot, "O31698/18-71", "SS", example_str)

        @test_throws AssertionError setannotresidue!(annot, "O31698/18-71",
                                                    String(rand('A':'Z',51)), example_str)

        @test_throws AssertionError setannotresidue!(annot, "O31698/18-71", "Feature Name",
                                                    example_str)

        @test ncolumns(annot) == 37

        setannotfile!(annot, "AC", "PF00571")
        setannotcolumn!(annot, "SS_cons", example_str)
        setannotsequence!(annot, "O31698/88-139", "OS", "Bacillus subtilis")

        @test getannotfile(annot, "AC") == "PF00571"
        @test getannotcolumn(annot, "SS_cons") == example_str
        @test getannotsequence(annot, "O31698/88-139", "OS") == "Bacillus subtilis"
        @test getannotresidue(annot, "O31698/18-71", "SS") == example_str

        @test getannotfile(annot, "An", "Default") == "Default"
        @test getannotcolumn(annot, "Other", "Default") == "Default"
        @test getannotsequence(annot, "O31698/1-88", "OS", "Default") == "Default"
        @test getannotresidue(annot, "O31698/1-88", "SS", "Default") == "Default"

        @test ncolumns(annot) == 37

        @test_throws DimensionMismatch  setannotresidue!(annot,
                                        "O31698/18-71", "AS", "__*__")

        @test_throws DimensionMismatch  setannotcolumn!(annot, "SS_cons",
                                        "---CCCCCHHHHHHHHHHHHHEEEEEEEEEEEEEEEEEEH---")
    end

    @testset "Copy, deepcopy and empty!" begin

        annot = Annotations()
        setannotresidue!(annot,"O31698/18-71","SS","CCCHHHHHHHHHHHHHHHEEEEEEEEEEEEEEEEHHH")
        setannotfile!(annot, "AC", "PF00571")
        setannotcolumn!(annot, "SS_cons", "CCCCCHHHHHHHHHHHHHEEEEEEEEEEEEEEEEEEH")
        setannotsequence!(annot, "O31698/88-139", "OS", "Bacillus subtilis")

        copy_annot = copy(annot)
        @test copy_annot == annot
        empty!(copy_annot)
        @test ncolumns(annot) == 37
        @test ncolumns(copy_annot) == -1

        @test !isempty(annot)
        @test isempty(copy_annot)

        deepcopy_annot = deepcopy(annot)
        @test deepcopy_annot == annot
        empty!(deepcopy_annot)
        @test ncolumns(annot) == 37
        @test ncolumns(deepcopy_annot) == -1

        @test !isempty(annot)
        @test isempty(deepcopy_annot)
    end

    @testset "Filter" begin

        @testset "Filter helpers" begin
            str_col = "abcd"
            str_map = "11,12,13,14"
            selector = [4,3,1]

            @test MSA._filter(str_col, selector) == "dca"
            @test MSA._filter_mapping(str_map, selector) == "14,13,11"
        end

        annot = Annotations()
        setannotresidue!(annot,"O31698/18-71","SS","CCCHHHHHHHHHHHHHHHEEEEEEEEEEEEEEEEHHH")
        setannotfile!(annot, "AC", "PF00571")
        setannotcolumn!(annot, "SS_cons", "CCCCCHHHHHHHHHHHHHEEEEEEEEEEEEEEEEEEH")
        setannotsequence!(annot, "O31698/88-139", "OS", "Bacillus subtilis")

        @test_throws AssertionError filtercolumns!(annot, [true, false, true])

        #filtersequences!(annot, IndexedArray(["O31698/88-139", "O31698/18-71"]), [false, true])
        #@test length( getannotsequence(annot) ) == 0
        #filtersequences!(annot, IndexedArray(["O31698/88-139", "O31698/18-71"]), [true, false])
        #@test length( getannotresidue(annot) ) == 0

        mask = collect("CCCCCHHHHHHHHHHHHHEEEEEEEEEEEEEEEEEEH") .!= Ref('E')
        filtercolumns!(annot, mask)
        @test ncolumns(annot) == 19
        @test getannotcolumn(annot, "SS_cons") == "CCCCCHHHHHHHHHHHHHH"

        filtercolumns!(annot, [1,2,19])
        @test ncolumns(annot) == 3
        @test getannotcolumn(annot, "SS_cons") == "CCH"
    end

    @testset "_rename_sequences" begin
        # Create an Annotations object with specific annotations
        annot = Annotations()
        setannotresidue!(annot, "O31698/18-71", "SS",  "CCCHHHHHHHHHHHHHHHEEEEEEEEEEEEEEEEHHH")
        setannotresidue!(annot, "O31698/72-140", "SS", "HHHHCCCCCEEEEEEEECCCCCCHHHHHHHHHHHHHH")
        setannotfile!(annot, "AC", "PF00571")
        setannotcolumn!(annot, "SS_cons", "CCCCCHHHHHHHHHHHHHEEEEEEEEEEEEEEEEEEH")
        setannotsequence!(annot, "O31698/88-139", "OS", "Bacillus subtilis")
        setannotsequence!(annot, "O31698/20-80", "OS", "Bacillus subtilis")
    
        # Define the old to new sequence name mapping
        old2new = Dict("O31698/18-71" => "NewSeq1", "O31698/88-139" => "NewSeq2")
    
        # Call the function with the annotations and mapping
        new_annotations = MSA._rename_sequences(annot, old2new)
    
        # Test cases
        # Check if the sequence names are correctly updated
        @test getannotresidue(new_annotations, "NewSeq1", "SS") == "CCCHHHHHHHHHHHHHHHEEEEEEEEEEEEEEEEHHH"
        @test getannotsequence(new_annotations, "NewSeq2", "OS") == "Bacillus subtilis"
    
        # Ensure that file and column annotations are unchanged
        @test getannotfile(new_annotations, "AC") == "PF00571"
        @test getannotcolumn(new_annotations, "SS_cons") == "CCCCCHHHHHHHHHHHHHEEEEEEEEEEEEEEEEEEH"
    
        # Check that the old sequence names are no longer present
        @test isempty(getannotresidue(new_annotations, "O31698/18-71", "SS", ""))
        @test isempty(getannotsequence(new_annotations, "O31698/88-139", "OS", ""))
    
        # Check that sequences not in the mapping remain unchanged
        @test getannotresidue(new_annotations, "O31698/72-140", "SS") == "HHHHCCCCCEEEEEEEECCCCCCHHHHHHHHHHHHHH"
        @test getannotsequence(new_annotations, "O31698/20-80", "OS") == "Bacillus subtilis"
    end    
    
    @testset "merge" begin
        @testset "different sources" begin
            original = Annotations(OrderedDict("key1" => "value1"), Dict(), Dict(), Dict())
            source1 = Annotations(OrderedDict("key2" => "value2"), Dict(), Dict(), Dict())
            source2 = Annotations(OrderedDict("key3" => "value3"), Dict(), Dict(), Dict())

            # Test merge without modifying the original
            merged = merge(original, source1, source2)
            @test length(merged.file) == 3
            @test merged.file["key1"] == "value1"
            @test merged.file["key2"] == "value2"
            @test merged.file["key3"] == "value3"

            # Original should not be altered by merge
            @test length(original.file) == 1
            @test original.file["key1"] == "value1"

            # Test merge!
            merge!(original, source1, source2)
            @test original == merged
        end

        @testset "overlapping keys" begin
            original = Annotations(OrderedDict("key" => "targetValue"), Dict(), Dict(), Dict())
            source1 = Annotations(OrderedDict("key" => "source1Value"), Dict(), Dict(), Dict())
            source2 = Annotations(OrderedDict("key" => "source2Value"), Dict(), Dict(), Dict())

            # Test merge without modifying the original
            merged = merge(original, source1, source2)
            @test length(merged.file) == 1
            @test merged.file["key"] == "source2Value"

            # Original should not be altered by merge
            @test length(original.file) == 1
            @test original.file["key"] == "targetValue"

            # Test merge!
            merge!(original, source1, source2)
            @test original == merged
        end

        @testset "empty Annotations" begin
            target = Annotations(OrderedDict("file_key1" => "file_value1"), Dict(), Dict(), Dict())
            empty_source = Annotations()

            merge!(target, empty_source)

            @test length(target.file) == 1
            @test target.file["file_key1"] == "file_value1"
            @test isempty(target.sequences)
            @test isempty(target.columns)
            @test isempty(target.residues)
        end
    
        @testset "partial overlap" begin
            target = Annotations(OrderedDict("key1" => "value1"), Dict(), Dict(), Dict())
            source = Annotations(OrderedDict("key1" => "new_value1", "key2" => "value2"), Dict(), Dict(), Dict())
    
            merge!(target, source)
    
            @test length(target.file) == 2
            @test target.file["key1"] == "new_value1"
            @test target.file["key2"] == "value2"
        end
    
        @testset "all annotation fields" begin
            target = Annotations(OrderedDict("file_key" => "file_value"), 
                                 Dict(("sequence_name", "annot_name") => "sequence_value"), 
                                 Dict("column_key" => "column_value"), 
                                 Dict(("residue_name", "residue_annot") => "residue_value"))
            source = Annotations(OrderedDict("file_key" => "new_file_value"), 
                                 Dict(("sequence_name", "annot_name") => "new_sequence_value"), 
                                 Dict("column_key" => "new_column_value"), 
                                 Dict(("residue_name", "residue_annot") => "new_residue_value"))
    
            merge!(target, source)
    
            @test target.file["file_key"] == "new_file_value"
            @test target.sequences[("sequence_name", "annot_name")] == "new_sequence_value"
            @test target.columns["column_key"] == "new_column_value"
            @test target.residues[("residue_name", "residue_annot")] == "new_residue_value"
        end
    end

    @testset "_rename_sequences" begin
        # Create a sample Annotations object with predefined annotations
        annotations = Annotations()
        setannotresidue!(annotations, "Seq1", "AnnotType1", "Value1")
        setannotresidue!(annotations, "Seq2", "AnnotType2", "Value2")
        setannotsequence!(annotations, "Seq1", "AnnotType3", "Value3")
        setannotsequence!(annotations, "Seq3", "AnnotType4", "Value4")
        setannotfile!(annotations, "FileKey", "FileValue")
        setannotcolumn!(annotations, "ColumnKey", "ValueX")
    
        # Define the old to new sequence name mapping
        old2new = Dict("Seq1" => "NewSeq1", "Seq2" => "NewSeq2")
    
        # Call the function with the annotations and mapping
        new_annotations = MSA._rename_sequences(annotations, old2new)
    
        # Test cases
        # Check if the sequence names are correctly updated
        @test getannotresidue(new_annotations, "NewSeq1", "AnnotType1") == "Value1"
        @test getannotresidue(new_annotations, "NewSeq2", "AnnotType2") == "Value2"
        @test getannotsequence(new_annotations, "NewSeq1", "AnnotType3") == "Value3"
    
        # Check if sequence names not in the mapping are retained
        @test getannotsequence(new_annotations, "Seq3", "AnnotType4") == "Value4"
    
        # Ensure that file and column annotations are unchanged
        @test getannotfile(new_annotations, "FileKey") == "FileValue"
        @test getannotcolumn(new_annotations, "ColumnKey") == "ValueX"
    end
end
