@testset "hcat" begin
    simple = joinpath(pwd(), "data", "simple.fasta")
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
end