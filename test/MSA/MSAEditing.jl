@testset "MSAEditing.jl" begin

    msa_types = (
        Matrix{Residue},
        NamedResidueMatrix{Array{Residue,2}},
        MultipleSequenceAlignment,
        AnnotatedMultipleSequenceAlignment
        )

    pf09645_sto = joinpath(pwd(), "data", "PF09645_full.stockholm")
    gaoetal2011 = joinpath(pwd(), "data", "Gaoetal2011.fasta")

    gaoetal_msas = [ read(gaoetal2011, FASTA, T) for T in msa_types ]
    pfam_msas    = [ read(pf09645_sto, Stockholm, T) for T in msa_types ]

    pfam    = pfam_msas[end]
    pfam_na = pfam_msas[end-1] # na: not annotated

    @testset "Boolean mask array" begin

        matrix_mask = pfam[1:1,:] .== GAP
        @test size(matrix_mask) == (1,110)

        @test size(MSA._column_mask(matrix_mask, pfam)) == (110,)
        @test MSA._column_mask(vec(matrix_mask), pfam) == vec(matrix_mask)
        @test MSA._column_mask(col -> col[1] == GAP, pfam) ==
              MSA._column_mask(matrix_mask, pfam)

        matrix_mask = pfam[:,1:1] .== Residue('-')
        @test size(matrix_mask) == (4,1)

        @test size(MSA._sequence_mask(matrix_mask, pfam)) == (4,)
        @test MSA._sequence_mask(vec(matrix_mask), pfam) == vec(matrix_mask)
        @test MSA._sequence_mask(seq -> seq[1] == GAP, pfam) ==
              MSA._sequence_mask(matrix_mask, pfam)
    end

    @testset "filtersequences!" begin

        for msa in pfam_msas[3:4]
            @test filtersequences!(copy(msa), 1:4 .< 3) == filtersequences(msa, 1:4 .< 3)
            @test nsequences(msa) == 4
        end

        for msa in pfam_msas
            @test getsequence(filtersequences(
                msa, [1,2,3,4] .== 3), 1) == getsequence(msa, 3)
            @test getsequence(filtersequences(
                msa, [1,2,3,4] .> 2), 2) == getsequence(msa, 4)

            @test getsequence(filtersequences(
                msa, Bool[false,false,true,true] ), 2) == getsequence(msa, 4)

            @test_throws AssertionError filtersequences(msa, [1,2,3,4,5,6,7,8,9,10] .> 2)
            @test_throws AssertionError filtersequences(msa, [1,2,3] .> 2)
        end
    end

    @testset "filtercolumns!" begin

        for msa in pfam_msas[3:4]
            @test filtercolumns!(copy(msa), 1:110 .< 11) == filtercolumns(msa, 1:110 .< 11)
            @test ncolumns(msa) == 110
        end

        for msa in pfam_msas
            @test vec(getresidues(getsequence(filtercolumns(
                msa, collect(1:110) .<= 10), 4))) == res"QTLNSYKMAE"

            @test_throws AssertionError filtercolumns(msa, [1,2,3] .> 2)
            @test_throws AssertionError filtercolumns(msa, collect(1:200) .<= 10)
        end

        @testset "filtercolumns! for sequences" begin

            annseq = getsequence(pfam, 4)
            seq = getsequence(pfam_na, 4) # na: not annotated

            # Sequences are matrices
            @test annseq[1, vec(annseq .== Residue('Q'))] == res"QQQQQQQQQQQQQQ"
            @test seq[1, vec(seq .== Residue('Q'))] == res"QQQQQQQQQQQQQQ"

            filtered_annseq = filtercolumns!(copy(annseq), annseq .== Residue('Q'))
            filtered_seq = filtercolumns!(copy(seq), seq .== Residue('Q'))

            @test filtercolumns(seq, seq .== Residue('Q')) ==
                filtercolumns(seq, seq .== Residue('Q'))

            @test_throws AssertionError filtercolumns(annseq, 1:(length(seq)-10) .> 2)
            @test_throws AssertionError filtercolumns(seq, 1:(length(seq)+10) .<= 10)

            # Sequences are matrices
            @test vec(getresidues(filtered_annseq)) == res"QQQQQQQQQQQQQQ"
            @test vec(getresidues(filtered_seq)) == res"QQQQQQQQQQQQQQ"

            @test getannotcolumn(filtered_annseq, "SS_cons") == "XHHEXXXXXXXXXX"
            @test getannotresidue(filtered_annseq,
                "F112_SSV1/3-112", "SS") == "XHHEXXXXXXXXXX"
        end
    end

    @testset "Reference" begin

        @testset "setreference" begin

            for msa in pfam_msas[2:4]
                copy_msa = copy(msa)
                setreference!(copy_msa, 4)
                # Sequences are matrices
                @test vec(getresidues(getsequence(copy_msa,1))) ==
                    res"QTLNSYKMAEIMYKILEKKGELTLEDILAQFEISVPSAYNIQRALKAICERHPDECEVQYKNRKTTFKWIKQEQKEEQKQEQTQDNIAKIFDAQPANFEQTDQGFIKAKQ"
                setreference!(copy_msa, "C3N734_SULIY/1-95")
                @test_throws ErrorException setreference!(copy_msa, "FALSE_SEQID/1-95")
                @test vec(getresidues(getsequence(copy_msa,4))) ==
                    res"QTLNSYKMAEIMYKILEKKGELTLEDILAQFEISVPSAYNIQRALKAICERHPDECEVQYKNRKTTFKWIKQEQKEEQKQEQTQDNIAKIFDAQPANFEQTDQGFIKAKQ"
                @test copy_msa == msa
            end
        end

        @testset "gapstrip!" begin

            for msa in pfam_msas[3:4]
                copy_msa = copy(msa)

                setreference!(copy_msa, 4)
                gapstrip!(copy_msa, gaplimit=1.0, coveragelimit=0.0)
                @test size(copy_msa) == (4, 110)

                setreference!(copy_msa, 1)
                gapstrip!(copy_msa)
                residuefraction(copy_msa[1,:]) == 1.0
            end

            for msa in pfam_msas[1:2]
                @test size(gapstrip(msa, gaplimit=1.0, coveragelimit=0.0)) == (4, 92)
                @test residuefraction(gapstrip(msa)[1,:]) == 1.0
            end
        end

        @testset "adjustreference!" begin

            for msa in pfam_msas[3:4]
                copy_msa = copy(msa)
                adjustreference!(copy_msa)
                residuefraction(copy_msa[1,:]) == 1.0
            end

            for msa in pfam_msas[1:2]
                @test residuefraction(adjustreference(msa)[1,:]) == 1.0
            end
        end
    end
end
