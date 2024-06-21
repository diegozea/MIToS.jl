@testset "hcat" begin
    simple = joinpath(DATA, "simple.fasta")
    msa = read_file(simple, FASTA, generatemapping=true)
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
            out_msa = read_file(tmp_file, Stockholm)
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
    msa = read_file(joinpath(DATA, "simple.fasta"), FASTA, generatemapping=true)
    msa2 = read_file(joinpath(DATA, "Gaoetal2011.fasta"), FASTA, generatemapping=true)

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

@testset "_find_gaps" begin
    @test MSA._find_gaps([2, 5, 6, 7, 8], 10) == [(2, 0), (5, 2), (11, 8)]
end

@testset "join MSAs" begin
    msa = read_file(joinpath(DATA, "simple.fasta"), FASTA, generatemapping=true)
    msa2 = read_file(joinpath(DATA, "Gaoetal2011.fasta"), FASTA, generatemapping=true)

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
            @test sum(startswith.(columnnames(gapped_msa), "gap:")) == 3 # the gap columns
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
        @test vec(gapped_col_seq) == Residue['A', '-', 'R']
    end

    @testset "_renumber_sequence_gaps" begin
        # Create a mock MSA
        M = rand(Residue, 5, 7)
        seqnames = ["seq1", "gap:3", "seq2", "gap:1", "gap:2"]
        column_names = ["col1", "col2", "col3", "col4", "col5", "col6", "col7"]
        named_matrix = MSA._namedresiduematrix(M, seqnames, column_names)
        annot = Annotations()
        setannotsequence!(annot, "gap:1", "AnnotationType1", "AnnotationValue1")
        setannotsequence!(annot, "seq1", "AnnotationType2", "AnnotationValue2")
        msa = AnnotatedMultipleSequenceAlignment(named_matrix, annot)
    
        # Apply the _renumber_sequence_gaps function
        new_msa = MSA._renumber_sequence_gaps(msa)
    
        # Test if the gap sequences are renumbered correctly
        @test sequencenames(new_msa) == ["seq1", "gap:1", "seq2", "gap:2", "gap:3"]
    
        # Verify that annotations are correctly transferred
        # "gap:1" in original becomes "gap:2" in new MSA
        @test getannotsequence(new_msa, "gap:2", "AnnotationType1") == "AnnotationValue1"
        # now, there is no annotations for "gap:1"
        @test isempty(getannotsequence(new_msa, "gap:1", "AnnotationType1", ""))
        # "seq1" remains unchanged
        @test getannotsequence(new_msa, "seq1", "AnnotationType2") == "AnnotationValue2"
    
        # Verify that the rest remains the same
        @test getresidues(new_msa) == getresidues(msa)
        @test columnnames(new_msa) == columnnames(msa)
    end

    @testset "_renumber_column_gaps" begin
        # Create a mock MSA
        M = rand(Residue, 2, 5)
        seqnames = ["seq1", "seq2"]
        column_names = ["col1", "gap:3", "col2", "gap:1", "gap:2"]
        named_matrix = MSA._namedresiduematrix(M, seqnames, column_names)
        msa = AnnotatedMultipleSequenceAlignment(named_matrix, Annotations())

        # Apply the _renumber_column_gaps function
        new_msa = MSA._renumber_column_gaps(msa)

        # Test if the gap columns are renumbered correctly
        @test columnnames(new_msa) == ["col1", "gap:1", "col2", "gap:2", "gap:3"]
    end

    @testset "_insert_sorted_gaps" begin 
        
        # NOTE: The default block_position is :before and the default axis is 1

        msa62 = msa2[:, 1:2]
        msa26 = msa2[1:2, :]
        
        @testset "gap sequences at the beginning and at the end" begin
            #
            #  1 -
            #  2 -
            # (3,1)
            # (4,2)
            # (5,3)
            # (6,4)
            #  - 5
            #  - 6
            #
            a = MSA._insert_sorted_gaps(msa62, msa62, [3, 4, 5, 6], [1, 2, 3, 4])
            @test size(a) == (8, 2)
            @test all(a[1:6, :] .!= GAP)
            @test all(a[7:8, :] .== GAP)
            b = MSA._insert_sorted_gaps(msa62, msa62, [1, 2, 3, 4], [3, 4, 5, 6])
            @test size(b) == (8, 2)
            @test all(b[1:2, :] .== GAP)
            @test all(b[3:8, :] .!= GAP)
        end

        @testset "gap columns at the beginning and at the end" begin
            a = MSA._insert_sorted_gaps(msa26, msa26, [3, 4, 5, 6], [1, 2, 3, 4], axis=2)
            @test size(a) == (2, 8)
            @test all(a[:, 1:6] .!= GAP)
            @test all(a[:, 7:8] .== GAP)
            b = MSA._insert_sorted_gaps(msa26, msa26, [1, 2, 3, 4], [3, 4, 5, 6], axis=2)
            @test size(b) == (2, 8)
            @test all(b[:, 1:2] .== GAP)
            @test all(b[:, 3:8] .!= GAP)
        end

        @testset "the unique matches are between the first and the last sequence" begin
            #
            # (1,1)
            #  2 -
            #  3 -
            #  4 -
            #  5 -
            #  - 2
            #  - 3
            #  - 4
            #  - 5
            # (6,6)
            #
            a = MSA._insert_sorted_gaps(msa62, msa62, [1, 6], [1, 6], block_position=:after)
            @test size(a) == (10, 2)
            @test all(a[1:5, :] .!= GAP)
            @test all(a[6:9, :] .== GAP)
            @test all(a[10:10, :] .!= GAP)
            b = MSA._insert_sorted_gaps(msa62, msa62, [1, 6], [1, 6])
            @test size(b) == (10, 2)
            @test all(b[1:1, :] .!= GAP)
            @test all(b[2:5, :] .== GAP)
            @test all(b[6:10, :] .!= GAP)
        end
        
        @testset "the unique matches are between the first and the last column" begin
            a = MSA._insert_sorted_gaps(msa26, msa26, [1, 6], [1, 6], axis=2, 
                block_position=:after)
            @test size(a) == (2, 10)
            @test all(a[:, 1:5] .!= GAP)
            @test all(a[:, 6:9] .== GAP)
            @test all(a[:, 10:10] .!= GAP)
            b = MSA._insert_sorted_gaps(msa26, msa26, [1, 6], [1, 6], axis=2) # default block_position=:before
            @test size(b) == (2, 10)
            @test all(b[:, 1:1] .!= GAP)
            @test all(b[:, 2:5] .== GAP)
            @test all(b[:, 6:10] .!= GAP)
        end

        @testset "gap sequences, gap sequences everywhere" begin
            #
            #  - 1
            # (1,2)
            #  - 3
            #  - 4
            # (2,5)
            #  3 -
            #  4 -
            #  5 -
            #  6 -
            #  - 6
            #
            a = MSA._insert_sorted_gaps(msa62, msa62, [1, 2], [2, 5], block_position=:after)
            @test size(a) == (10, 2)
            @test all(a[1:1, :] .== GAP)
            @test all(a[2:2, :] .!= GAP)
            @test all(a[3:4, :] .== GAP)
            @test all(a[5:9, :] .!= GAP)
            @test all(a[10:10, :] .== GAP)
            b = MSA._insert_sorted_gaps(msa62, msa62, [2, 5], [1, 2])
            @test size(b) == (10, 2)
            @test all(b[1:5, :] .!= GAP)
            @test all(b[6:9, :] .== GAP)
            @test all(b[10:10, :] .!= GAP)
        end

        @testset "gap columns, gap columns everywhere" begin
            a = MSA._insert_sorted_gaps(msa26, msa26, [1, 2], [2, 5], axis=2, 
                block_position=:after)
            @test size(a) == (2, 10)
            @test all(a[:, 1:1] .== GAP)
            @test all(a[:, 2:2] .!= GAP)
            @test all(a[:, 3:4] .== GAP)
            @test all(a[:, 5:9] .!= GAP)
            @test all(a[:, 10:10] .== GAP)
            b = MSA._insert_sorted_gaps(msa26, msa26, [2, 5], [1, 2], axis=2)
            @test size(b) == (2, 10)
            @test all(b[:, 1:5] .!= GAP)
            @test all(b[:, 6:9] .== GAP)
            @test all(b[:, 10:10] .!= GAP)
        end
    end

    @testset "_add_gaps_in_b" begin
        msa62 = msa2[:, 1:2]  # msa62 for sequences test
        msa26 = msa2[1:2, :]  # msa26 for columns test
    
        @testset "gaps at the beginning and at the end" begin
            # 1 2 3 4 5 6 - -
            # - - 1 2 3 4 5 6

            # a: 1 2 3 4 5 6 
            # b: - - 1 2 3 4 

            # sequences
            b = MSA._add_gaps_in_b(msa62, msa62, 3:6, 1:4)
            @test size(b) == (6, 2)
            @test all(b[1:2, :] .== GAP)
            @test all(b[3:6, :] .!= GAP)
            # columns
            b = MSA._add_gaps_in_b(msa26, msa26, 3:6, 1:4, 2)
            @test size(b) == (2, 6)
            @test all(b[:, 1:2] .== GAP)
            @test all(b[:, 3:6] .!= GAP)

            # a: 3 4 5 6 - -
            # b: 1 2 3 4 5 6

            # sequences
            a = MSA._add_gaps_in_b(msa62, msa62, 1:4, 3:6)
            @test size(a) == (6, 2)
            @test all(a[1:4, :] .!= GAP)
            @test all(a[5:6, :] .== GAP)

            # columns
            a = MSA._add_gaps_in_b(msa26, msa26, 1:4, 3:6, 2)
            @test size(a) == (2, 6)
            @test all(a[:, 1:4] .!= GAP)
            @test all(a[:, 5:6] .== GAP)
        end

        @testset "unique matches between first and last sequence/column" begin
            # 1 2 3 4 5 - - - - 6
            # 1 - - - - 2 3 4 5 6
    
            # a: 1 2 3 4 5 6
            # b: 1 - - - - 6

            # sequences
            b_seq = MSA._add_gaps_in_b(msa62, msa62, [1, 6], [1, 6])
            @test size(b_seq) == (6, 2)
            @test all(b_seq[1, :] .!= GAP)
            @test all(b_seq[2:5, :] .== GAP)
            @test all(b_seq[6, :] .!= GAP)
    
            # columns
            b_col = MSA._add_gaps_in_b(msa26, msa26, [1, 6], [1, 6], 2)
            @test size(b_col) == (2, 6)
            @test all(b_col[:, 1] .!= GAP)
            @test all(b_col[:, 2:5] .== GAP)
            @test all(b_col[:, 6] .!= GAP)
        end

        @testset "gap sequences, gap sequences everywhere" begin
            # - 1 - - 2 3 4 5 6 -
            # 1 2 3 4 5 - - - - 6

            # a: 1 2 3 4 5 6
            # b: 2 5 - - - -
            # sequences
            b_seq = MSA._add_gaps_in_b(msa62, msa62, [1, 2], [2, 5])
            @test size(b_seq) == (6, 2)
            @test all(b_seq[1:2, :] .!= GAP)
            @test all(b_seq[3:6, :] .== GAP)

            # columns
            b_col = MSA._add_gaps_in_b(msa26, msa26, [1, 2], [2, 5], 2)
            @test size(b_col) == (2, 6)
            @test all(b_col[:, 1:2] .!= GAP)
            @test all(b_col[:, 3:6] .== GAP)

            # a: - 1 - - 2 -
            # b: 1 2 3 4 5 6

            # sequences
            a_seq = MSA._add_gaps_in_b(msa62, msa62, [2, 5], [1, 2])
            @test size(a_seq) == (6, 2)
            @test all(a_seq[1, :] .== GAP)
            @test all(a_seq[2, :] .!= GAP)
            @test all(a_seq[3:4, :] .== GAP)
            @test all(a_seq[5, :] .!= GAP)
            @test all(a_seq[6, :] .== GAP)

            # columns
            a_col = MSA._add_gaps_in_b(msa26, msa26, [2, 5], [1, 2], 2)
            @test size(a_col) == (2, 6)
            @test all(a_col[:, 1] .== GAP)
            @test all(a_col[:, 2] .!= GAP)
            @test all(a_col[:, 3:4] .== GAP)
            @test all(a_col[:, 5] .!= GAP)
            @test all(a_col[:, 6] .== GAP)
        end
    end

    @testset "join MSAs" begin
        msa62 = msa2[:, 1:2]  # msa62 for sequences test
        msa26 = msa2[1:2, :]  # msa26 for columns test

        @testset "gaps at the beginning and at the end" begin
            # a: 1 2 3 4 5 6 - -
            # b: - - 1 2 3 4 5 6

            @testset "inner join" begin
                # a: 3 4 5 6
                # b: 1 2 3 4

                @testset "sequences" begin
                    ab = join(msa62, msa62, 3:6 .=> 1:4, kind=:inner, axis=1)
                    @test size(ab) == (4, 4)
                    @test all(ab .!= GAP)
                    @test sequencenames(ab) == ["SEQ3_&_SEQ1", "SEQ4_&_SEQ2", "SEQ5_&_SEQ3", 
                        "SEQ6_&_SEQ4"]
                    @test columnnames(ab) == ["1_1", "1_2", "2_1", "2_2"]
                end

                @testset "columns" begin 
                    ab = join(msa26, msa26, 3:6 .=> 1:4, kind=:inner, axis=2)
                    @test size(ab) == (4, 4)
                    @test all(ab .!= GAP)
                    @test sequencenames(ab) == ["1_SEQ1", "1_SEQ2", "2_SEQ1", "2_SEQ2"]
                    @test columnnames(ab) == ["3", "4", "5", "6"]
                end
            end

            @testset "outer join" begin
                # a: 1 2 3 4 5 6 - -
                # b: - - 1 2 3 4 5 6

                @testset "sequences" begin
                    ab = join(msa62, msa62, 3:6 .=> 1:4) #  default: kind=:outer, axis=1
                    @test size(ab) == (8, 4)
                    @test sum(ab .== GAP, dims=1) == [2 2 2 2]
                    @test vec(sum(ab .== GAP, dims=2)) == [2, 2, 0, 0, 0, 0, 2, 2]
                    @test sequencenames(ab) == ["SEQ1_&_gap:1", "SEQ2_&_gap:2", 
                        "SEQ3_&_SEQ1", "SEQ4_&_SEQ2", "SEQ5_&_SEQ3", "SEQ6_&_SEQ4", 
                        "gap:1_&_SEQ5", "gap:2_&_SEQ6"]
                    @test columnnames(ab) == ["1_1", "1_2", "2_1", "2_2"]
                end

                @testset "columns" begin
                    ab = join(msa26, msa26, 3:6 .=> 1:4, axis=2) # default: kind=:outer
                    @test size(ab) == (4, 8)
                    @test vec(sum(ab .== GAP, dims=2)) == [2, 2, 2, 2]
                    @test vec(sum(ab .== GAP, dims=1)) == [2, 2, 0, 0, 0, 0, 2, 2]
                    @test sequencenames(ab) == ["1_SEQ1", "1_SEQ2", "2_SEQ1", "2_SEQ2"]
                    @test columnnames(ab) == ["1", "2", "3", "4", "5", "6", "gap:1", "gap:2"]
                end
            end

            @testset "left join" begin
                # a: 1 2 3 4 5 6
                # b: - - 1 2 3 4

                @testset "sequences" begin
                    ab = join(msa62, msa62, 3:6 .=> 1:4, kind=:left, axis=1)
                    @test size(ab) == (6, 4)
                    @test sum(ab.==GAP, dims=1) == [0 0 2 2]
                    @test vec(sum(ab.==GAP, dims=2)) == [2, 2, 0, 0, 0, 0]
                    @test sequencenames(ab) == ["SEQ1_&_gap:1", "SEQ2_&_gap:2", 
                        "SEQ3_&_SEQ1", "SEQ4_&_SEQ2", "SEQ5_&_SEQ3", "SEQ6_&_SEQ4"]
                    @test columnnames(ab) == ["1_1", "1_2", "2_1", "2_2"]
                end

                @testset "columns" begin
                    ab = join(msa26, msa26, 3:6 .=> 1:4, kind=:left, axis=2)
                    @test size(ab) == (4, 6)
                    @test vec(sum(ab.==GAP, dims=2)) == [0, 0, 2, 2]
                    @test vec(sum(ab.==GAP, dims=1)) == [2, 2, 0, 0, 0, 0]
                    @test sequencenames(ab) == ["1_SEQ1", "1_SEQ2", "2_SEQ1", "2_SEQ2"]
                    @test columnnames(ab) == ["1", "2", "3", "4", "5", "6"]
                end
            end

            @testset "right join" begin
                # a: 3 4 5 6 - -
                # b: 1 2 3 4 5 6

                @testset "sequences" begin
                    ab = join(msa62, msa62, 3:6 .=> 1:4, kind=:right, axis=1)
                    @test size(ab) == (6, 4)
                    @test sum(ab.==GAP, dims=1) == [2 2 0 0]
                    @test vec(sum(ab.==GAP, dims=2)) == [0, 0, 0, 0, 2, 2]
                    @test sequencenames(ab) == ["SEQ3_&_SEQ1", "SEQ4_&_SEQ2", 
                        "SEQ5_&_SEQ3", "SEQ6_&_SEQ4", "gap:1_&_SEQ5", "gap:2_&_SEQ6"]
                end

                @testset "columns" begin
                    ab = join(msa26, msa26, 3:6 .=> 1:4, kind=:right, axis=2)
                    @test size(ab) == (4, 6)
                    @test vec(sum(ab.==GAP, dims=2)) == [2, 2, 0, 0]
                    @test vec(sum(ab.==GAP, dims=1)) == [0, 0, 0, 0, 2, 2]
                    @test sequencenames(ab) == ["1_SEQ1", "1_SEQ2", "2_SEQ1", "2_SEQ2"]
                    @test columnnames(ab) == ["3", "4", "5", "6", "gap:1", "gap:2"]
                end
            end
        end

        @testset "unique matches between first and last sequence/column" begin
            # 1 2 3 4 5 - - - - 6
            # 1 - - - - 2 3 4 5 6

            @testset "inner join" begin
                # a: 1 6
                # b: 1 6

                @testset "sequences" begin
                    ab = join(msa62, msa62, [1, 6] .=> [1, 6], kind=:inner, axis=1)
                    @test size(ab) == (2, 4)
                    @test all(ab .!= GAP)
                    @test sequencenames(ab) == ["SEQ1", "SEQ6"]
                    @test columnnames(ab) == ["1_1", "1_2", "2_1", "2_2"]
                    @test ab == Residue['D' 'A' 'D' 'A'; 'D' 'A' 'D' 'A']
                end

                @testset "columns" begin
                    ab = join(msa26, msa26, [1, 6] .=> [1, 6], kind=:inner, axis=2)
                    @test size(ab) == (4, 2)
                    @test all(ab .!= GAP)
                    @test sequencenames(ab) == ["1_SEQ1", "1_SEQ2", "2_SEQ1", "2_SEQ2"]
                    @test columnnames(ab) == ["1", "6"]
                    @test ab == Residue['D' 'E'; 'D' 'F'; 'D' 'E'; 'D' 'F']
                end
            end

            @testset "left join" begin
                # a: 1 2 3 4 5 6
                # b: 1 - - - - 6

                @testset "sequences" begin
                    ab = join(msa62, msa62, [1, 6] .=> [1, 6], kind=:left)
                    @test size(ab) == (6, 4)
                    @test sum(ab.==GAP, dims=1) == [0 0 4 4]
                    @test vec(sum(ab.==GAP, dims=2)) == [0, 2, 2, 2, 2, 0]
                    @test sequencenames(ab) == ["SEQ1", "SEQ2", "SEQ3", 
                        "SEQ4", "SEQ5", "SEQ6"]
                    @test columnnames(ab) == ["1_1", "1_2", "2_1", "2_2"]
                end

                @testset "columns" begin
                    ab = join(msa26, msa26, [1, 6] .=> [1, 6], kind=:left, axis=2)
                    @test size(ab) == (4, 6)
                    @test vec(sum(ab.==GAP, dims=2)) == [0, 0, 4, 4]
                    @test vec(sum(ab.==GAP, dims=1)) == [0, 2, 2, 2, 2, 0]
                    @test sequencenames(ab) == ["1_SEQ1", "1_SEQ2", "2_SEQ1", "2_SEQ2"]
                    @test columnnames(ab) == ["1", "2", "3", "4", "5", "6"]
                end
            end

            @testset "right join" begin
                # a: 1 - - - - 6
                # b: 1 2 3 4 5 6

                @testset "sequences" begin
                    ab = join(msa62, msa62, [1, 6] .=> [1, 6], kind=:right)
                    @test size(ab) == (6, 4)
                    @test sum(ab.==GAP, dims=1) == [4 4 0 0]
                    @test vec(sum(ab.==GAP, dims=2)) == [0, 2, 2, 2, 2, 0]
                    @test sequencenames(ab) == ["SEQ1", "SEQ2", "SEQ3", 
                        "SEQ4", "SEQ5", "SEQ6"]
                    @test columnnames(ab) == ["1_1", "1_2", "2_1", "2_2"]
                end

                @testset "columns" begin
                    ab = join(msa26, msa26, [1, 6] .=> [1, 6], kind=:right, axis=2)
                    @test size(ab) == (4, 6)
                    @test vec(sum(ab.==GAP, dims=2)) == [4, 4, 0, 0]
                    @test vec(sum(ab.==GAP, dims=1)) == [0, 2, 2, 2, 2, 0]
                    @test sequencenames(ab) == ["1_SEQ1", "1_SEQ2", "2_SEQ1", "2_SEQ2"]
                    @test columnnames(ab) == ["1", "gap:1", "gap:2", "gap:3", "gap:4", "6"]
                end
            end

            @testset "outer join" begin
                # a: 1 2 3 4 5 - - - - 6
                # b: 1 - - - - 2 3 4 5 6

                @testset "sequences" begin
                    ab = join(msa62, msa62, [(1, 1), (6, 6)], kind=:outer, axis=1)
                    @test size(ab) == (10, 4)
                    @test sum(ab.==GAP, dims=1) == [4 4 4 4]
                    @test vec(sum(ab.==GAP, dims=2)) == [0, 2, 2, 2, 2, 2, 2, 2, 2, 0]

                    @test sequencenames(ab) == ["SEQ1", "SEQ2_&_gap:1", "SEQ3_&_gap:2", 
                        "SEQ4_&_gap:3", "SEQ5_&_gap:4", "gap:1_&_SEQ2", "gap:2_&_SEQ3",
                        "gap:3_&_SEQ4", "gap:4_&_SEQ5", "SEQ6"]
                    @test columnnames(ab) == ["1_1", "1_2", "2_1", "2_2"]

                    # kind=:outer and axis=1 are the default values
                    @test ab == join(msa62, msa62, [(1, 1), (6, 6)])

                    @testset "pairing types: _find_pairing_positions" begin
                        # This is mostly to test the _find_pairing_positions function
                        # that it is used by join at the beginning of the function.

                        # Positions
                        # ---------
                        # Vector of pairs
                        @test ab == join(msa62, msa62, [1 => 1, 6 => 6])
                        # Vector of tuples
                        @test ab == join(msa62, msa62, [(1, 1), (6, 6)])
                        # Vector of vectors
                        @test ab == join(msa62, msa62, [[1, 1], [6, 6]])
                        # Vector of named tuples
                        @test ab == join(msa62, msa62, [(a=1, b=1), (a=6, b=6)])
                        # Tuple of pairs
                        @test ab == join(msa62, msa62, (1 => 1, 6 => 6))
                        # Tuple of tuples
                        @test ab == join(msa62, msa62, ((1, 1), (6, 6)))
                        # Tuple of vectors
                        @test ab == join(msa62, msa62, ([1, 1], [6, 6]))
                        # Tuple of named tuples
                        @test ab == join(msa62, msa62, ((a=1, b=1), (a=6, b=6)))
                        
                        # Sequence names
                        # --------------
                        # Vector of pairs
                        @test ab == join(msa62, msa62, ["SEQ1" => "SEQ1", "SEQ6" => "SEQ6"])
                        # Vector of tuples
                        @test ab == join(msa62, msa62, [("SEQ1", "SEQ1"), ("SEQ6", "SEQ6")])
                        # Vector of vectors
                        @test ab == join(msa62, msa62, [["SEQ1", "SEQ1"], ["SEQ6", "SEQ6"]])
                        # Vector of named tuples
                        @test ab == join(msa62, msa62, [(a="SEQ1", b="SEQ1"), (a="SEQ6", b="SEQ6")])
                        # Tuple of pairs
                        @test ab == join(msa62, msa62, ("SEQ1" => "SEQ1", "SEQ6" => "SEQ6"))
                        # Tuple of tuples
                        @test ab == join(msa62, msa62, (("SEQ1", "SEQ1"), ("SEQ6", "SEQ6")))
                        # Tuple of vectors
                        @test ab == join(msa62, msa62, (["SEQ1", "SEQ1"], ["SEQ6", "SEQ6"]))
                        # Tuple of named tuples
                        @test ab == join(msa62, msa62, ((a="SEQ1", b="SEQ1"), (a="SEQ6", b="SEQ6")))
                        
                        # OrderedDict
                        # -----------
                        # NOTE: join change its behavior depending on whether the pairing
                        # positions are sorted or not. Here, OrderedDict ensures that the 
                        # positions are sorted.
                        @test ab == join(msa62, msa62, OrderedDict{Int, Int}(1 => 1, 6 => 6))
                    end
                end

                @testset "columns" begin
                    ab = join(msa26, msa26, [(1, 1), (6, 6)], kind=:outer, axis=2)
                    @test size(ab) == (4, 10)
                    @test vec(sum(ab.==GAP, dims=2)) == [4, 4, 4, 4]
                    @test sum(ab.==GAP, dims=1) == [0 2 2 2 2 2 2 2 2 0]
                    @test sequencenames(ab) == ["1_SEQ1", "1_SEQ2", "2_SEQ1", "2_SEQ2"]
                    @test columnnames(ab) == ["1", "2", "3", "4", "5", "gap:1", "gap:2", 
                        "gap:3", "gap:4", "6"]
                end
            end
        end

        @testset "gap sequences, gap sequences everywhere" begin
            # - 1 - - 2 3 4 5 6 -
            # 1 2 3 4 5 - - - - 6

            @testset "inner join" begin
                # a: 1 2
                # b: 2 5

                @testset "sequences" begin
                    ab = join(msa62, msa62, [1, 2] .=> [2, 5], kind=:inner, axis=1)
                    @test size(ab) == (2, 4)
                    @test all(ab .!= GAP)
                    @test sequencenames(ab) == ["SEQ1_&_SEQ2", "SEQ2_&_SEQ5"]
                    @test columnnames(ab) == ["1_1", "1_2", "2_1", "2_2"]
                    @test ab == Residue['D' 'A' 'D' 'A'; 'D' 'A' 'D' 'A']
                    @testset "annotations" begin
                        @test getannotfile(ab, "NCol") == "6_&_6"
                        @test gethcatmapping(ab) == [1, 1, 2, 2]
                        for seq in 1:2
                            @test getsequencemapping(ab, seq) == [1, 2, 1, 2]
                        end
                    end
                end

                @testset "columns" begin
                    ab = join(msa26, msa26, [1, 2] .=> [2, 5], kind=:inner, axis=2)
                    @test size(ab) == (4, 2)
                    @test all(ab .!= GAP)
                    @test sequencenames(ab) == ["1_SEQ1", "1_SEQ2", "2_SEQ1", "2_SEQ2"]
                    @test columnnames(ab) == ["1", "2"]
                    @test ab == Residue['D' 'A'; 'D' 'A'; 'A' 'E'; 'A' 'E']
                    @testset "annotations" begin
                        @test getannotfile(ab, "1_NCol") == "6"
                        @test getannotfile(ab, "2_NCol") == "6"
                        @test getannotfile(ab, "1_ColMap") == "1,2"
                        @test getannotfile(ab, "2_ColMap") == "2,5"
                        @test_throws ErrorException gethcatmapping(ab)    
                        @test getsequencemapping(ab, 1) == [1, 2]
                        @test getsequencemapping(ab, 3) == [2, 5]
                    end
                end
            end

            @testset "left join" begin
                # a: 1 2 3 4 5 6
                # b: 2 5 - - - -

                @testset "sequences" begin
                    ab = join(msa62, msa62, [1, 2] .=> [2, 5], kind=:left)
                    @test size(ab) == (6, 4)
                    @test sum(ab.==GAP, dims=1) == [0 0 4 4]
                    @test vec(sum(ab.==GAP, dims=2)) == [0, 0, 2, 2, 2, 2]
                    @test sequencenames(ab) == ["SEQ1_&_SEQ2", "SEQ2_&_SEQ5", "SEQ3", 
                        "SEQ4", "SEQ5_&_gap:1", "SEQ6"]
                    @test columnnames(ab) == ["1_1", "1_2", "2_1", "2_2"]
                    @testset "annotations" begin
                        @test getannotfile(ab, "NCol") == "6_&_6"
                        @test gethcatmapping(ab) == [1, 1, 2, 2]
                        @test getannotfile(ab, "ColMap") == "1,2,1,2"
                        for seq in 1:2
                            @test getsequencemapping(ab, seq) == [1, 2, 1, 2]
                        end
                        for seq in 3:6
                            @test getsequencemapping(ab, seq) == [1, 2, 0, 0]
                        end
                    end
                end

                @testset "columns" begin
                    ab = join(msa26, msa26, [1, 2] .=> [2, 5], kind=:left, axis=2)
                    @test size(ab) == (4, 6)
                    @test vec(sum(ab.==GAP, dims=2)) == [0, 0, 4, 4]
                    @test vec(sum(ab.==GAP, dims=1)) == [0, 0, 2, 2, 2, 2]
                    @test sequencenames(ab) == ["1_SEQ1", "1_SEQ2", "2_SEQ1", "2_SEQ2"]
                    @test columnnames(ab) == ["1", "2", "3", "4", "5", "6"]
                    @testset "annotations" begin
                        @test getannotfile(ab, "1_NCol") == "6"
                        @test getannotfile(ab, "2_NCol") == "6"
                        @test getannotfile(ab, "1_ColMap") == "1,2,3,4,5,6"
                        @test getannotfile(ab, "2_ColMap") == "2,5,,,,"
                        @test_throws ErrorException gethcatmapping(ab)
                        for seq in 1:2
                            @test getsequencemapping(ab, seq) == [1, 2, 3, 4, 5, 6]
                        end
                        for seq in 3:4
                            @test getsequencemapping(ab, seq) == [2, 5, 0, 0, 0, 0]
                        end
                    end
                end
            end

            @testset "right join" begin
                # a: - 1 - - 2 -
                # b: 1 2 3 4 5 6

                @testset "sequences" begin
                    ab = join(msa62, msa62, [1, 2] .=> [2, 5], kind=:right)
                    @test size(ab) == (6, 4)
                    @test sum(ab.==GAP, dims=1) == [4 4 0 0]
                    @test vec(sum(ab.==GAP, dims=2)) == [2, 0, 2, 2, 0, 2]
                    @test sequencenames(ab) == ["gap:1_&_SEQ1", "SEQ1_&_SEQ2", "SEQ3", 
                        "SEQ4", "SEQ2_&_SEQ5", "SEQ6"]
                end

                @testset "columns" begin
                    ab = join(msa26, msa26, [1, 2] .=> [2, 5], kind=:right, axis=2)
                    @test size(ab) == (4, 6)
                    @test vec(sum(ab.==GAP, dims=2)) == [4, 4, 0, 0]
                    @test vec(sum(ab.==GAP, dims=1)) == [2, 0, 2, 2, 0, 2]
                    @test sequencenames(ab) == ["1_SEQ1", "1_SEQ2", "2_SEQ1", "2_SEQ2"]
                    @test columnnames(ab) == ["gap:1", "1", "gap:2", "gap:3", "2", "gap:4"]
                end
            end

            @testset "outer join" begin
                # a: - 1 - - 2 3 4 5 6 -
                # b: 1 2 3 4 5 - - - - 6

                @testset "sequences" begin
                    ab = join(msa62, msa62, [1, 2] .=> [2, 5], kind=:outer)
                    @test size(ab) == (10, 4)
                    @test sum(ab.==GAP, dims=1) == [4 4 4 4]
                    @test vec(sum(ab.==GAP, dims=2)) == [2, 0, 2, 2, 0, 2, 2, 2, 2, 2]
                    @test sequencenames(ab) == ["gap:1_&_SEQ1", "SEQ1_&_SEQ2", "gap:2_&_SEQ3",
                        "gap:3_&_SEQ4", "SEQ2_&_SEQ5", "SEQ3_&_gap:1", "SEQ4_&_gap:2",
                        "SEQ5_&_gap:3", "SEQ6_&_gap:4", "gap:4_&_SEQ6"]
                    @test columnnames(ab) == ["1_1", "1_2", "2_1", "2_2"]

                    @testset "pairing types: _find_pairing_positions" begin
                        # See the notes in the previous _find_pairing_positions testset.

                        # Positions
                        # ---------
                        # Vector of pairs
                        @test ab == join(msa62, msa62, [1 => 2, 2 => 5])
                        # Vector of tuples
                        @test ab == join(msa62, msa62, [(1, 2), (2, 5)])
                        # Vector of vectors
                        @test ab == join(msa62, msa62, [[1, 2], [2, 5]])
                        # Vector of named tuples
                        @test ab == join(msa62, msa62, [(a=1, b=2), (a=2, b=5)])
                        # Tuple of pairs
                        @test ab == join(msa62, msa62, (1 => 2, 2 => 5))
                        # Tuple of tuples
                        @test ab == join(msa62, msa62, ((1, 2), (2, 5)))
                        # Tuple of vectors
                        @test ab == join(msa62, msa62, ([1, 2], [2, 5]))
                        # Tuple of named tuples
                        @test ab == join(msa62, msa62, ((a=1, b=2), (a=2, b=5)))
                        
                        # Sequence names
                        # --------------
                        # Vector of pairs
                        @test ab == join(msa62, msa62, ["SEQ1" => "SEQ2", "SEQ2" => "SEQ5"])
                        # Vector of tuples
                        @test ab == join(msa62, msa62, [("SEQ1", "SEQ2"), ("SEQ2", "SEQ5")])
                        # Vector of vectors
                        @test ab == join(msa62, msa62, [["SEQ1", "SEQ2"], ["SEQ2", "SEQ5"]])
                        # Vector of named tuples
                        @test ab == join(msa62, msa62, [(a="SEQ1", b="SEQ2"), (a="SEQ2", b="SEQ5")])
                        # Tuple of pairs
                        @test ab == join(msa62, msa62, ("SEQ1" => "SEQ2", "SEQ2" => "SEQ5"))
                        # Tuple of tuples
                        @test ab == join(msa62, msa62, (("SEQ1", "SEQ2"), ("SEQ2", "SEQ5")))
                        # Tuple of vectors
                        @test ab == join(msa62, msa62, (["SEQ1", "SEQ2"], ["SEQ2", "SEQ5"]))
                        # Tuple of named tuples
                        @test ab == join(msa62, msa62, ((a="SEQ1", b="SEQ2"), (a="SEQ2", b="SEQ5")))

                        # OrderedDict
                        # -----------
                        @test ab == join(msa62, msa62, OrderedDict{Int, Int}(1 => 2, 2 => 5))
                        @test ab == join(msa62, msa62, OrderedDict{String, String}("SEQ1" => "SEQ2", "SEQ2" => "SEQ5"))
                    end

                    @testset "annotations" begin
                        @test getannotfile(ab, "NCol") == "6_&_6"
                        @test gethcatmapping(ab) == [1, 1, 2, 2]
                        @test getannotfile(ab, "ColMap") == "1,2,1,2"
                        for seq in [2, 5]
                            @test getsequencemapping(ab, seq) == [1, 2, 1, 2]
                        end
                        for seq in [1, 3, 4, 10]
                            @test getsequencemapping(ab, seq) == [0, 0, 1, 2]
                        end
                        for seq in [6, 7, 8, 9]
                            @test getsequencemapping(ab, seq) == [1, 2, 0, 0]
                        end
                    end
                end

                @testset "columns" begin
                    ab = join(msa26, msa26, [1, 2] .=> [2, 5], kind=:outer, axis=2)
                    @test size(ab) == (4, 10)
                    @test vec(sum(ab.==GAP, dims=2)) == [4, 4, 4, 4]
                    @test vec(sum(ab.==GAP, dims=1)) == [2, 0, 2, 2, 0, 2, 2, 2, 2, 2]
                    @test sequencenames(ab) == ["1_SEQ1", "1_SEQ2", "2_SEQ1", "2_SEQ2"]
                    @test columnnames(ab) == ["gap:1", "1", "gap:2", "gap:3", "2", "3", 
                        "4", "5", "6", "gap:4"]
                    @testset "annotations" begin
                        @test getannotfile(ab, "1_NCol") == "6"
                        @test getannotfile(ab, "2_NCol") == "6"
                        @test getannotfile(ab, "1_ColMap") == ",1,,,2,3,4,5,6,"
                        @test getannotfile(ab, "2_ColMap") == "1,2,3,4,5,,,,,6"
                        @test_throws ErrorException gethcatmapping(ab)
                        @test getsequencemapping(ab, 1) == [0, 1, 0, 0, 2, 3, 4, 5, 6, 0]
                        @test getsequencemapping(ab, 3) == [1, 2, 3, 4, 5, 0, 0, 0, 0, 6]
                    end
                end
            end
        end

        @testset "unsorted positions" begin
            # When the pairing positions are not sorted, the join function will behave 
            # differently depending on the kind of join. For an inner join,
            # the join function will keep the order indicating by the pairing.
            # For the outer join it will first match the positions as indicated in the 
            # pairing input. Then, it will add the remaining unmatched sequences/columns 
            # for `msa_a` and `msa_b` in that order. In the case of left and right join,
            # the order of the sequences/columns will be the same as the order in the 
            # left/right MSA.

        
            # - 1 4 2 3 5 6 - -
            # 1 3 4 2 - - - 5 6

            @testset "inner join" begin
                # a: 1 4 2
                # b: 3 4 2

                @testset "sequences" begin
                    ab = join(msa62, msa62, [1, 4, 2] .=> [3, 4, 2], kind=:inner)
                    @test size(ab) == (3, 4)
                    @test all(ab .!= GAP)
                    @test sequencenames(ab) == ["SEQ1_&_SEQ3", "SEQ4", "SEQ2"]
                    @test columnnames(ab) == ["1_1", "1_2", "2_1", "2_2"]
                end

                @testset "columns" begin
                    ab = join(msa26, msa26, [1, 4, 2] .=> [3, 4, 2], kind=:inner, axis=2)
                    @test size(ab) == (4, 3)
                    @test all(ab .!= GAP)
                    @test sequencenames(ab) == ["1_SEQ1", "1_SEQ2", "2_SEQ1", "2_SEQ2"]
                    @test columnnames(ab) == ["1", "4", "2"]
                end 
            end

            @testset "left join" begin
                # a: 1 2 3 4 5 6
                # b: 3 2 - 4 - -

                @testset "sequences" begin
                    ab = join(msa62, msa62, [1, 4, 2] .=> [3, 4, 2], kind=:left)
                    @test size(ab) == (6, 4)
                    @test sum(ab.==GAP, dims=1) == [0 0 3 3]
                    @test vec(sum(ab.==GAP, dims=2)) == [0, 0, 2, 0, 2, 2]
                    @test sequencenames(ab) == ["SEQ1_&_SEQ3", "SEQ2", 
                        "SEQ3_&_gap:1", "SEQ4", "SEQ5", "SEQ6"]
                    @test columnnames(ab) == ["1_1", "1_2", "2_1", "2_2"]
                end

                @testset "columns" begin
                    ab = join(msa26, msa26, [1, 4, 2] .=> [3, 4, 2], kind=:left, axis=2)
                    @test size(ab) == (4, 6)
                    @test vec(sum(ab.==GAP, dims=2)) == [0, 0, 3, 3]
                    @test vec(sum(ab.==GAP, dims=1)) == [0, 0, 2, 0, 2, 2]
                    @test sequencenames(ab) == ["1_SEQ1", "1_SEQ2", "2_SEQ1", "2_SEQ2"]
                    @test columnnames(ab) == ["1", "2", "3", "4", "5", "6"]
                end
            end

            @testset "right join" begin
                # a: - 2 1 4 - -
                # a: 1 2 3 4 5 6

                @testset "sequences" begin
                    ab = join(msa62, msa62, [1, 4, 2] .=> [3, 4, 2], kind=:right)
                    @test size(ab) == (6, 4)
                    @test sum(ab.==GAP, dims=1) == [3 3 0 0]
                    @test vec(sum(ab.==GAP, dims=2)) == [2, 0, 0, 0, 2, 2]
                    @test sequencenames(ab) == ["gap:1_&_SEQ1", "SEQ2", "SEQ1_&_SEQ3", 
                        "SEQ4", "SEQ5", "SEQ6"]
                end

                @testset "columns" begin
                    ab = join(msa26, msa26, [1, 4, 2] .=> [3, 4, 2], kind=:right, axis=2)
                    @test size(ab) == (4, 6)
                    @test vec(sum(ab.==GAP, dims=2)) == [3, 3, 0, 0]
                    @test vec(sum(ab.==GAP, dims=1)) == [2, 0, 0, 0, 2, 2]
                    @test sequencenames(ab) == ["1_SEQ1", "1_SEQ2", "2_SEQ1", "2_SEQ2"]
                    @test columnnames(ab) == ["gap:1", "2", "1", "4", "gap:2", "gap:3"]
                end
            end

            @testset "outer join" begin
                # a: 1 4 2 3 5 6 - - - 
                # b: 3 4 2 - - - 1 5 6

                @testset "sequences" begin
                    ab = join(msa62, msa62, [1, 4, 2] .=> [3, 4, 2], kind=:outer)
                    @test size(ab) == (9, 4)
                    @test sum(ab.==GAP, dims=1) == [3 3 3 3]
                    @test vec(sum(ab.==GAP, dims=2)) == [0, 0, 0, 2, 2, 2, 2, 2, 2]
                    @test sequencenames(ab) == ["SEQ1_&_SEQ3", "SEQ4", "SEQ2", 
                        "SEQ3_&_gap:1", "SEQ5_&_gap:2", "SEQ6_&_gap:3", "gap:1_&_SEQ1",
                        "gap:2_&_SEQ5", "gap:3_&_SEQ6"]
                    @test columnnames(ab) == ["1_1", "1_2", "2_1", "2_2"]
                    @testset "annotations" begin
                        @test getannotfile(ab, "NCol") == "6_&_6"
                        @test gethcatmapping(ab) == [1, 1, 2, 2]
                        @test getannotfile(ab, "ColMap") == "1,2,1,2"
                        for seq in 1:3
                            @test getsequencemapping(ab, seq) == [1, 2, 1, 2]
                        end
                        for seq in 4:6
                            @test getsequencemapping(ab, seq) == [1, 2, 0, 0]
                        end
                        for seq in 7:9
                            @test getsequencemapping(ab, seq) == [0, 0, 1, 2]
                        end
                    end
                end

                @testset "columns" begin
                    ab = join(msa26, msa26, [1, 4, 2] .=> [3, 4, 2], kind=:outer, axis=2)
                    @test size(ab) == (4, 9)
                    @test vec(sum(ab.==GAP, dims=2)) == [3, 3, 3, 3]
                    @test vec(sum(ab.==GAP, dims=1)) == [0, 0, 0, 2, 2, 2, 2, 2, 2]
                    @test sequencenames(ab) == ["1_SEQ1", "1_SEQ2", "2_SEQ1", "2_SEQ2"]
                    @test columnnames(ab) == ["1", "4", "2", "3", "5", "6", "gap:1", 
                        "gap:2", "gap:3"]
                    @testset "annotations" begin
                        @test getannotfile(ab, "1_NCol") == "6"
                        @test getannotfile(ab, "2_NCol") == "6"
                        @test getannotfile(ab, "1_ColMap") == "1,4,2,3,5,6,,,"
                        @test getannotfile(ab, "2_ColMap") == "3,4,2,,,,1,5,6"
                        @test_throws ErrorException gethcatmapping(ab)
                        for seq in 1:2
                            @test getsequencemapping(ab, seq) == [1, 4, 2, 3, 5, 6, 0, 0, 0]
                        end
                        for seq in 3:4
                            @test getsequencemapping(ab, seq) == [3, 4, 2, 0, 0, 0, 1, 5, 6]
                        end
                    end
                end
            end
        end

        @testset "ArgumentErrors" begin
            # axis is not 1 or 2
            @test_throws ArgumentError join(msa62, msa62, [1, 2] .=> [3, 4], axis=3)
            # kind is not :inner, :left, :right, or :outer
            @test_throws ArgumentError join(msa62, msa62, [1, 2] .=> [3, 4], kind=:iner)
            # pairing is empty
            @test_throws ArgumentError join(msa62, msa62, [])
            # each element of the pairing is not a pair
            @test_throws ArgumentError join(msa62, msa62, [1, 2, 3])
            # the list of positions are not of the same length
            @test_throws ArgumentError join(msa62, msa62, [1, 2], [3, 4, 5])
        end

        @testset "Using two position lists instead of a list of pairs" begin
            for axis in [1, 2]
                for kind in [:inner, :left, :right, :outer]
                    @test join(msa62, msa62, [1, 2], [2, 1], kind=kind, axis=axis) ==
                        join(msa62, msa62, [(1, 2), (2, 1)], kind=kind, axis=axis)
                end
            end
        end
    end
end
