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

@testset "fuse" begin
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
end