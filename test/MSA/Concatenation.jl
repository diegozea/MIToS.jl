@testset "hcat" begin
    simple = joinpath(DATA, "simple.fasta")
    msa = read(simple, FASTA, generatemapping=true)
    setannotresidue!(msa, "ONE", "example", "ab")
    setannotresidue!(msa, "TWO", "example", "cd")
    setannotresidue!(msa, "ONE", "OnlyONE", "xx")
    setannotresidue!(msa, "TWO", "OnlyTWO", "yy")
    msa_2 = copy(msa)
    setannotcolumn!(msa_2, "example", "HE")

    same_ref = Residue[
        'A' 'R' 'A' 'R'
        'R' 'A' 'R' 'A'
    ]

    diff_ref = Residue[
        'A' 'R' 'R' 'A'
        'R' 'A' 'A' 'R'
    ]

    annot_same = hcat(msa, msa_2)
    annot_diff = hcat(msa_2, msa[[2, 1], :])
    msa_same = hcat(
        MultipleSequenceAlignment(msa), 
        MultipleSequenceAlignment(msa))
    msa_diff = hcat(
        MultipleSequenceAlignment(msa), 
        MultipleSequenceAlignment(msa)[[2, 1], :])

    @testset "Same sequence names" begin
        for concatenated_msa in (annot_same, msa_same)
            @test size(concatenated_msa) == (2, 4)
            @test concatenated_msa == same_ref
            @test sequencenames(concatenated_msa) == ["ONE", "TWO"]
            @test columnnames(concatenated_msa) == ["1_1", "1_2", "2_1", "2_2"]
            @test getcolumnmapping(concatenated_msa) == [1, 2, 1, 2]
            if concatenated_msa isa AnnotatedMultipleSequenceAlignment
                @test getsequencemapping(concatenated_msa, "ONE") == [1, 2, 1, 2]
                @test getsequencemapping(concatenated_msa, "TWO") == [1, 2, 1, 2]
                @test getannotresidue(concatenated_msa, "ONE", "example") == "abab"
                @test getannotresidue(concatenated_msa, "TWO", "example") == "cdcd"
                @test getannotresidue(concatenated_msa, "ONE", "OnlyONE") == "xxxx"
                @test getannotresidue(concatenated_msa, "TWO", "OnlyTWO") == "yyyy"
                @test getannotcolumn(concatenated_msa, "example") == "  HE"
            end
            @test gethcatmapping(concatenated_msa) == [1, 1, 2, 2]
        end
    end

    @testset "Different sequence names" begin
        for concatenated_msa in (annot_diff, msa_diff)
            @test size(concatenated_msa) == (2, 4)
            @test concatenated_msa == diff_ref
            @test sequencenames(concatenated_msa) == ["ONE_&_TWO", "TWO_&_ONE"]
            @test columnnames(concatenated_msa) == ["1_1", "1_2", "2_1", "2_2"]
            @test getcolumnmapping(concatenated_msa) == [1, 2, 1, 2]
            if concatenated_msa isa AnnotatedMultipleSequenceAlignment
                @test getsequencemapping(concatenated_msa, "ONE_&_TWO") == [1, 2, 1, 2]
                @test getsequencemapping(concatenated_msa, "TWO_&_ONE") == [1, 2, 1, 2]
                @test getannotresidue(concatenated_msa, "ONE_&_TWO", "example") == "abcd"
                @test getannotresidue(concatenated_msa, "TWO_&_ONE", "example") == "cdab"
                @test getannotresidue(concatenated_msa, "ONE_&_TWO", "OnlyONE") == "xx  "
                @test getannotresidue(concatenated_msa, "ONE_&_TWO", "OnlyTWO") == "  yy"
                @test getannotresidue(concatenated_msa, "TWO_&_ONE", "OnlyONE") == "  xx"
                @test getannotresidue(concatenated_msa, "TWO_&_ONE", "OnlyTWO") == "yy  "
                @test getannotcolumn(concatenated_msa, "example") == "HE  "
            end
            @test gethcatmapping(concatenated_msa) == [1, 1, 2, 2]
        end
    end

    @testset "Multiple concatenations and annot column/residue" begin
        concatenated_msa = hcat(msa, msa_2, msa, msa_2, msa)
        @test getannotcolumn(concatenated_msa, "example") == "  HE  HE  "
    end

    @testset "IO" begin
        path = tempdir()
        tmp_file = joinpath(path, ".tmp.stockholm")
        try
            write(tmp_file, annot_diff, Stockholm)
            out_msa = read(tmp_file, Stockholm)
            @test columnnames(out_msa) == columnnames(annot_diff)
            @test out_msa == annot_diff
        finally
            if isfile(tmp_file)
                rm(tmp_file)
            end
        end
    end

    @testset "Inception" begin
        concatenated_in = hcat(msa, msa_2)
        concatenated_diff_a = hcat(msa[[2, 1], :], msa_2)
        concatenated_diff_b = hcat(msa_2, msa[[2, 1], :])

        @testset "concatenated concatenated" begin
            concatenated_out = hcat(concatenated_in, concatenated_in)
            concat_ab = hcat(concatenated_diff_a, concatenated_diff_b)

            @test size(concatenated_out) == (2, 8)
            @test sequencenames(concatenated_out) == ["ONE", "TWO"]
            @test columnnames(concatenated_out) == [
                "1_1", "1_2", "2_1", "2_2", "3_1", "3_2", "4_1", "4_2"
                ]
            @test getcolumnmapping(concatenated_out) == [1, 2, 1, 2, 1, 2, 1, 2]
            @test getsequencemapping(concatenated_out, "ONE") == [1, 2, 1, 2, 1, 2, 1, 2]
            @test getsequencemapping(concatenated_out, "TWO") == [1, 2, 1, 2, 1, 2, 1, 2]
            @test getsequencemapping(concatenated_out, "ONE") == [1, 2, 1, 2, 1, 2, 1, 2]
            @test getsequencemapping(concatenated_out, "TWO") == [1, 2, 1, 2, 1, 2, 1, 2]
            @test getannotresidue(concatenated_out, "ONE", "example") == "abababab"
            @test getannotresidue(concatenated_out, "TWO", "example") == "cdcdcdcd"
            @test getannotresidue(concatenated_out, "ONE", "OnlyONE") == "xxxxxxxx"
            @test getannotresidue(concatenated_out, "TWO", "OnlyTWO") == "yyyyyyyy"
            @test getannotcolumn(concatenated_out, "example") == "  HE  HE"
            @test gethcatmapping(concatenated_out) == [1, 1, 2, 2, 3, 3, 4, 4]

            @test size(concat_ab) == (2, 8)
            @test sequencenames(concat_ab) == [
                "TWO_&_ONE_&_ONE_&_TWO", "ONE_&_TWO_&_TWO_&_ONE"]
            @test columnnames(concat_ab) == [
                "1_1", "1_2", "2_1", "2_2", "3_1", "3_2", "4_1", "4_2"]
            @test getcolumnmapping(concat_ab) == [1, 2, 1, 2, 1, 2, 1, 2]
            @test getsequencemapping(concat_ab, 
                "TWO_&_ONE_&_ONE_&_TWO") == [1, 2, 1, 2, 1, 2, 1, 2]
            @test getsequencemapping(concat_ab, 
                "ONE_&_TWO_&_TWO_&_ONE") == [1, 2, 1, 2, 1, 2, 1, 2]
            @test getannotresidue(concat_ab, 
                "TWO_&_ONE_&_ONE_&_TWO", "example") == "cdababcd"
            @test getannotresidue(concat_ab, 
                "ONE_&_TWO_&_TWO_&_ONE", "example") == "abcdcdab"
            @test getannotresidue(concat_ab, 
                "TWO_&_ONE_&_ONE_&_TWO", "OnlyONE") == "  xxxx  "
            @test getannotresidue(concat_ab, 
                "TWO_&_ONE_&_ONE_&_TWO", "OnlyTWO") == "yy    yy"
            @test getannotcolumn(concat_ab, "example") == "  HEHE  "
            @test gethcatmapping(concat_ab) == [1, 1, 2, 2, 3, 3, 4, 4]
        end

        @testset "concatenated non_concatenated" begin
            concatenated_msas = [
                hcat(concatenated_in, msa),
                hcat(msa, concatenated_in)
            ]
            for (i, concatenated_out) in enumerate(concatenated_msas)
                @test size(concatenated_out) == (2, 6)
                @test sequencenames(concatenated_out) == ["ONE", "TWO"]
                @test columnnames(concatenated_out) == ["1_1", "1_2", "2_1", "2_2", "3_1", "3_2"]
                @test getcolumnmapping(concatenated_out) == [1, 2, 1, 2, 1, 2]
                @test getsequencemapping(concatenated_out, "ONE") == [1, 2, 1, 2, 1, 2]
                @test getsequencemapping(concatenated_out, "TWO") == [1, 2, 1, 2, 1, 2]
                @test getsequencemapping(concatenated_out, "ONE") == [1, 2, 1, 2, 1, 2]
                @test getsequencemapping(concatenated_out, "TWO") == [1, 2, 1, 2, 1, 2]
                @test getannotresidue(concatenated_out, "ONE", "example") == "ababab"
                @test getannotresidue(concatenated_out, "TWO", "example") == "cdcdcd"
                @test getannotresidue(concatenated_out, "ONE", "OnlyONE") == "xxxxxx"
                @test getannotresidue(concatenated_out, "TWO", "OnlyTWO") == "yyyyyy"
                if i == 1
                    @test getannotcolumn(concatenated_out, "example") == "  HE  "
                else
                    @test getannotcolumn(concatenated_out, "example") == "    HE"
                end
                @test gethcatmapping(concatenated_out) == [1, 1, 2, 2, 3, 3]
            end

            concat_a = hcat(concatenated_diff_a, msa)

            @test size(concat_a) == (2, 6)
            @test sequencenames(concat_a) == [
                "TWO_&_ONE_&_ONE", "ONE_&_TWO_&_TWO"]
            @test columnnames(concat_a) == [
                "1_1", "1_2", "2_1", "2_2", "3_1", "3_2"]
            @test getcolumnmapping(concat_a) == [1, 2, 1, 2, 1, 2]
            @test getsequencemapping(concat_a, "TWO_&_ONE_&_ONE") == [1, 2, 1, 2, 1, 2]
            @test getsequencemapping(concat_a, "ONE_&_TWO_&_TWO") == [1, 2, 1, 2, 1, 2]
            @test getannotresidue(concat_a, "TWO_&_ONE_&_ONE", "example") == "cdabab"
            @test getannotresidue(concat_a, "ONE_&_TWO_&_TWO", "example") == "abcdcd"
            @test getannotresidue(concat_a, 
                "TWO_&_ONE_&_ONE", "OnlyONE") == "  xxxx"
            @test getannotresidue(concat_a, 
                "TWO_&_ONE_&_ONE", "OnlyTWO") == "yy    "
            @test getannotcolumn(concat_a, "example") == "  HE  "
            @test gethcatmapping(concat_a) == [1, 1, 2, 2, 3, 3]
        end
    end
end

@testset "vcat" begin
    msa = read(joinpath(DATA, "simple.fasta"), FASTA, generatemapping=true)
    msa2 = read(joinpath(DATA, "Gaoetal2011.fasta"), FASTA, generatemapping=true)

    @testset "seqnames" begin
        seqnames = sequencenames(msa)
        seqnames2 = sequencenames(msa2)

        new_names, label_map = MIToS.MSA._v_concatenated_seq_names(msa, msa)
        @test new_names == ["1_ONE", "1_TWO", "2_ONE", "2_TWO"]
        @test isempty(label_map) # the are no previous prefixes/labels
        new_names, label_map = MIToS.MSA._v_concatenated_seq_names(msa, msa, msa)
        @test new_names == ["1_ONE", "1_TWO", "2_ONE", "2_TWO", "3_ONE", "3_TWO"]
        @test isempty(label_map)
        new_names, label_map = MIToS.MSA._v_concatenated_seq_names(msa, msa2)
        @test new_names == vcat(["1_$n" for n in seqnames], ["2_$n" for n in seqnames2])
        @test isempty(label_map)
        new_names, label_map = MIToS.MSA._v_concatenated_seq_names(msa2, msa)
        @test new_names == vcat(["1_$n" for n in seqnames2], ["2_$n" for n in seqnames])
        @test isempty(label_map)
        
        @testset "Sequence name mapping" begin
            new_names, label_mapping = MIToS.MSA._v_concatenated_seq_names(msa, msa)
            mapping = MIToS.MSA._get_seqname_mapping_vcat(new_names, msa, msa)
            @test isempty(label_mapping) # the are no previous prefixes/labels
            @test length(mapping) == 4
            @test mapping[(1, "ONE")] == "1_ONE"
            @test mapping[(1, "TWO")] == "1_TWO"
            @test mapping[(2, "ONE")] == "2_ONE"
            @test mapping[(2, "TWO")] == "2_TWO"
        end
    end

    @testset "vcat examples" begin
        setannotresidue!(msa, "ONE", "res_example", "ab")
        setannotcolumn!(msa, "example", "HE")
        concatenated_11 = vcat(msa, msa)
        setannotfile!(msa, "file_example", "file annotation")
        setannotsequence!(msa, "ONE", "seq_example", "seq annotation")
        msa_2 = copy(msa)
        setannotfile!(msa_2, "file_example_msa2", "file annotation msa2")
        setannotsequence!(msa_2, "ONE", "seq_example_msa2", "seq annotation msa2")
        setannotresidue!(msa_2, "ONE", "res_example_msa2", "AB")
        setannotcolumn!(msa_2, "example_msa2", "he")
        concatenated = vcat(msa, msa_2)

        # matrix
        @test size(concatenated) == (4, 2)
        @test concatenated == Residue[
            'A' 'R'
            'R' 'A'
            'A' 'R'
            'R' 'A'
        ]

        # names
        @test sequencenames(concatenated) == ["1_ONE", "1_TWO", "2_ONE", "2_TWO"]
        @test columnnames(concatenated) == ["1", "2"]

        @testset "unannotated aligned objects" begin
           msa_unannot = MultipleSequenceAlignment(msa)
           vcat_unannot = vcat(msa_unannot, msa_unannot)
           @test size(vcat_unannot) == (4, 2)
            @test vcat_unannot == Residue[
                'A' 'R'
                'R' 'A'
                'A' 'R'
                'R' 'A'
            ]
            @test sequencenames(vcat_unannot) == ["1_ONE", "1_TWO", "2_ONE", "2_TWO"]
            @test columnnames(vcat_unannot) == ["1", "2"]
        end

        @testset "vcat annotations" begin
            # sequence annotations: disambiguated by msa number (prefix)
            @test getannotsequence(concatenated, "1_ONE", "seq_example") == "seq annotation"
            @test getannotsequence(concatenated, "2_ONE", "seq_example_msa2") == "seq annotation msa2"
            @test getsequencemapping(concatenated, "1_ONE") == [1, 2]
            @test getsequencemapping(concatenated, "2_ONE") == [1, 2]
            # residue annotations: disambiguated by sequence name
            @test getannotresidue(concatenated, "1_ONE", "res_example") == "ab"
            @test getannotresidue(concatenated, "2_ONE", "res_example") == "ab"
            @test getannotresidue(concatenated, "2_ONE", "res_example_msa2") == "AB"
            # column annotations: disambiguated by msa number (prefix)
            @test getannotcolumn(concatenated, "1_example") == "HE"
            @test getannotcolumn(concatenated, "2_example") == "HE"
            @test getannotcolumn(concatenated, "2_example_msa2") == "he"
        end

        @testset "Inception" begin
            concatenated_111 = vcat(msa, concatenated_11)

            @test size(concatenated_111) == (6, 2)
            @test concatenated_111 == Residue[
                'A' 'R'
                'R' 'A'
                'A' 'R'
                'R' 'A'
                'A' 'R'
                'R' 'A'
            ]
            @test sequencenames(concatenated_111) == [
                "1_ONE", "1_TWO", "2_ONE", "2_TWO", "3_ONE", "3_TWO"]
            @test columnnames(concatenated_111) == ["1", "2"]
            # sequence annotations: disambiguated by msa number (prefix)
            @test getannotsequence(concatenated_111, "1_ONE", "SeqMap") == "1,2"
            @test getannotsequence(concatenated_111, "2_ONE", "SeqMap") == "1,2"
            @test getannotsequence(concatenated_111, "3_ONE", "SeqMap") == "1,2"
            # residue annotations: disambiguated by sequence name
            @test getannotresidue(concatenated_111, "1_ONE", "res_example") == "ab"
            @test getannotresidue(concatenated_111, "2_ONE", "res_example") == "ab"
            @test getannotresidue(concatenated_111, "3_ONE", "res_example") == "ab"
            # file annotations: disambiguated by msa number (prefix)
            @test getannotfile(concatenated_111, "1_NCol") == "2"
            @test getannotfile(concatenated_111, "2_NCol") == "2"
            @test getannotfile(concatenated_111, "3_NCol") == "2"
            # column annotations: disambiguated by msa number (prefix)
            @test getannotcolumn(concatenated_111, "1_example") == "HE"
            @test getannotcolumn(concatenated_111, "2_example") == "HE"
            @test getannotcolumn(concatenated_111, "3_example") == "HE"
        end
    end
end

@testset "join MSAs" begin
    msa = read(joinpath(DATA, "simple.fasta"), FASTA, generatemapping=true)
    msa2 = read(joinpath(DATA, "Gaoetal2011.fasta"), FASTA, generatemapping=true)

    @testset "column gaps" begin
        h_gaps = MIToS.MSA._gap_columns(msa, 3)
        @test all(==(GAP), h_gaps)
        @test size(h_gaps) == (2, 3)
        h_concatenated = hcat(msa[:, 1:1], h_gaps, msa[:, 2:2])
        @test size(h_concatenated) == (2, 5)
        @test h_concatenated == Residue[
            'A' '-' '-' '-' 'R'
            'R' '-' '-' '-' 'A'
        ]
        @test sequencenames(h_concatenated) == ["ONE", "TWO"]
        @test getcolumnmapping(h_concatenated) == [1, 0, 0, 0, 2]
    end

    @testset "sequence gaps" begin
        v_gaps = MIToS.MSA._gap_sequences(msa, ["SEQ1", "SEQ2", "SEQ3"])
        @test all(==(GAP), v_gaps)
        @test size(v_gaps) == (3, 2)
        v_concatenated = vcat(msa[1:1, :], v_gaps, msa[2:2, :])
        @test size(v_concatenated) == (5, 2)
        @test v_concatenated == Residue[
            'A' 'R'
            '-' '-'
            '-' '-'
            '-' '-'
            'R' 'A'
        ]
        vcat_seqnames = sequencenames(v_concatenated)
        @test vcat_seqnames == ["1_ONE", "2_SEQ1", "2_SEQ2", "2_SEQ3", "3_TWO"]
        @test getsequencemapping(v_concatenated, "1_ONE") == [1, 2]
        @test getsequencemapping(v_concatenated, "2_SEQ1") == [0, 0]
    end

    @testset "insert gap sequences" begin 
        at_the_start = MIToS.MSA._insert_gap_sequences(msa, ["SEQ1", "SEQ2", "SEQ3"], 1)
        in_the_middle = MIToS.MSA._insert_gap_sequences(msa, ["SEQ1", "SEQ2", "SEQ3"], 2)
        at_the_end = MIToS.MSA._insert_gap_sequences(msa, ["SEQ1", "SEQ2", "SEQ3"], 3)
        
        for gapped_msa in [at_the_start, in_the_middle, at_the_end]
            @test size(gapped_msa) == (5, 2)
            @test sum(gapped_msa .== GAP) == 6 # 3 x 2
            @test sort(sequencenames(gapped_msa)) == ["ONE", "SEQ1", "SEQ2", "SEQ3", "TWO"]
        end

        # check that the annotations are the same
        delete_annotated_modifications!(at_the_start)
        delete_annotated_modifications!(in_the_middle)
        delete_annotated_modifications!(at_the_end)
        @test annotations(at_the_start) == annotations(at_the_end)
        @test annotations(at_the_start) == annotations(in_the_middle)

        # check the position of the gap blocks
        @test sum(at_the_start[1:3, :] .== GAP) == 6
        @test sum(in_the_middle[2:4, :] .== GAP) == 6
        @test sum(at_the_end[3:5, :] .== GAP) == 6
    end

    @testset "insert gap columns" begin
        at_the_start = MIToS.MSA._insert_gap_columns(msa, 3, 1)
        in_the_middle = MIToS.MSA._insert_gap_columns(msa, 3, 2)
        at_the_end = MIToS.MSA._insert_gap_columns(msa, 3, 3)

        for gapped_msa in [at_the_start, in_the_middle, at_the_end]
            @test size(gapped_msa) == (2, 5) # (2 , 2 + 3)
            @test sum(gapped_msa .== GAP) == 6 # 2 x 3
            @test sort(columnnames(gapped_msa))[1:2] == ["1", "2"] # the original columns
            @test sum(startswith.(columnnames(gapped_msa), "padding:")) == 3 # the gap columns
        end

        # Check the position of the gap columns
        @test sum(at_the_start[:, 1:3] .== GAP) == 6 # gap columns at start
        @test sum(in_the_middle[:, 2:4] .== GAP) == 6 # gap columns in the middle
        @test sum(at_the_end[:, 3:5] .== GAP) == 6 # gap columns at end

        # Check if original columns are preserved correctly
        @test at_the_start[:, [4, 5]] == msa
        @test in_the_middle[:, [1, 5]] == msa
        @test at_the_end[:, [1, 2]] == msa

        # Check for correct column mapping after inserting gaps
        @test getcolumnmapping(at_the_start) == [0, 0, 0, 1, 2]
        @test getcolumnmapping(in_the_middle) == [1, 0, 0, 0, 2]
        @test getcolumnmapping(at_the_end) == [1, 2, 0, 0, 0]

        # Check for correct sequence mapping after inserting gaps
        @test getsequencemapping(at_the_start, "ONE") == [0, 0, 0, 1, 2]
        @test getsequencemapping(in_the_middle, "ONE") == [1, 0, 0, 0, 2]
        @test getsequencemapping(at_the_end, "ONE") == [1, 2, 0, 0, 0]

        # Check annotations are preserved correctly; this operation should not add or 
        # delete annotations
        @test length(annotations(at_the_start)) == length(annotations(msa))
        @test length(annotations(in_the_middle)) == length(annotations(msa))
        @test length(annotations(at_the_end)) == length(annotations(msa))

        @testset "MSA without NCol annotation" begin
            msa_no_NCol = deepcopy(msa)
            annotfile = getannotfile(msa_no_NCol)
            delete!(annotfile, "NCol")
            @test getannotfile(msa_no_NCol, "NCol", "") == ""
    
            # test that no NCol annotation is added
            gapped_msa = MIToS.MSA._insert_gap_columns(msa_no_NCol, 3, 1)
            @test !haskey(annotations(gapped_msa).file, "NCol")
        end
    end

    @testset "Insert gaps: Unannotated MSAs" begin
        ref_gapped_msa_seq = MIToS.MSA._insert_gap_sequences(msa, ["SEQ1", "SEQ2"], 3)
        ref_gapped_msa_col = MIToS.MSA._insert_gap_columns(msa, 3, 3)

        for msa_unannot in [MultipleSequenceAlignment(msa), namedmatrix(msa)]
            gapped_msa_seq = MIToS.MSA._insert_gap_sequences(msa_unannot, ["SEQ1", "SEQ2"], 3)
            gapped_msa_col = MIToS.MSA._insert_gap_columns(msa_unannot, 3, 3)

            # tests to ensure sequences and columns are inserted correctly
            @test gapped_msa_seq == ref_gapped_msa_seq
            @test gapped_msa_col == ref_gapped_msa_col

            # test for equal sequence and column names
            @test sequencenames(gapped_msa_seq) == sequencenames(ref_gapped_msa_seq)
            @test columnnames(gapped_msa_seq) == columnnames(ref_gapped_msa_seq)
            @test sequencenames(gapped_msa_col) == sequencenames(ref_gapped_msa_col)
            @test columnnames(gapped_msa_col)[1:2] == columnnames(ref_gapped_msa_col)[1:2]
        end
    end

    @testset "Insert gaps: Invalid gap positions" begin
        # Position less than 1
        @test_throws ArgumentError MIToS.MSA._insert_gap_sequences(msa, ["SEQ1", "SEQ2"], 0)
        @test_throws ArgumentError MIToS.MSA._insert_gap_columns(msa, 3, 0)

        # Position greater than the number of sequences/columns are valid and used to 
        # insert gaps at the end of the MSA
    end

    @testset "Insert gaps: Single sequence/column MSAs" begin
        single_seq_msa = msa[1:1, :]
        single_col_msa = msa[:, 1:1]

        # Test for a single sequence
        gapped_seq_seq = MIToS.MSA._insert_gap_sequences(single_seq_msa, ["SEQ1"], 2) # at the end
        @test size(gapped_seq_seq) == (2, ncolumns(msa))
        @test gapped_seq_seq == Residue['A' 'R'; '-' '-']
        
        gapped_seq_col = MIToS.MSA._insert_gap_columns(single_seq_msa, 1, 2) # at the middle
        @test size(gapped_seq_col) == (1, 3)
        @test gapped_seq_col == Residue[
                'A' '-' 'R'
            ]

        # Test for a single column
        gapped_col_col = MIToS.MSA._insert_gap_columns(single_col_msa, 1, 2) # at the end
        @test size(gapped_col_col) == (nsequences(msa), 2)
        @test gapped_col_col == Residue['A' '-'; 'R' '-']

        gapped_col_seq = MIToS.MSA._insert_gap_sequences(single_col_msa, ["SEQ1"], 2) # at the middle
        @test size(gapped_col_seq) == (3, 1)
        @test gapped_col_seq == Residue['A'; '-'; 'R';;]
    end

    @testset "_compress_array!" begin
        @test MIToS.MSA._compress_array!([1, 2, 3, 6, 7, 8, 10, 20, 21, 22]) == [
            1:3, 6:8, 10:10, 20:22]
        @test MIToS.MSA._compress_array!(Int[]) == UnitRange{Int}[]
        @test MIToS.MSA._compress_array!([5]) == [5:5]
        @test MIToS.MSA._compress_array!([2, 1, 3, 4, 5]) == [2:2, 1:1, 3:5]
        @test MIToS.MSA._compress_array!([1, 3, 5, 7, 9]) == [1:1, 3:3, 5:5, 7:7, 9:9]
        @test MIToS.MSA._compress_array!([4, 4, 4, 4]) == [4:4]
    end

end
