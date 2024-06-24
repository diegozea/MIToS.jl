@testset "Shuffle" begin

    msa_types = (
        Matrix{Residue},
        NamedResidueMatrix{Array{Residue,2}},
        MultipleSequenceAlignment,
        AnnotatedMultipleSequenceAlignment
        )

    N = length(msa_types)

    pf09645_sto = joinpath(DATA, "PF09645_full.stockholm")

    msas = [ read_file(pf09645_sto, Stockholm, T) for T in msa_types ]
    gaps = [ msa .== GAP for msa in msas ]
    lcol = [ mean(msa .== Residue('L'), dims=1) for msa in msas ]
    lseq = [ mean(msa .== Residue('L'), dims=2) for msa in msas ]

    Random.seed!(42)

    @testset "General" begin

        msa = msas[1]

        for dim in [1,2]

            aln = copy(msa)
            shuffle_msa!(aln, dims=dim)
            @test aln != msa
            @test (aln .== GAP) == gaps[1] # default: fixed gaps

            aln = copy(msa)
            shuffle_msa!(aln, dims=dim, fixedgaps=true)
            @test aln != msa
            @test (aln .== GAP) == gaps[1]

            aln = copy(msa)
            shuffle_msa!(aln, dims=dim, fixedgaps=false)
            @test aln != msa
            @test (aln .== GAP) != gaps[1]
        end

        @test_throws AssertionError shuffle_msa(msa, dims=0, fixedgaps=true)
        @test_throws AssertionError shuffle_msa(msa, dims=3, fixedgaps=true)
    end

    @testset "Columns" begin

        for i in 1:N
            # Fixed gaps
            msa = msas[i]
            aln = shuffle_msa(msa, dims=2, fixedgaps=true)
            @test aln != getresidues(msa)
            @test (aln .== GAP) == gaps[i]
            @test lcol[i] == mean(aln .== Residue('L'), dims=1)
            @test lseq[i] != mean(aln .== Residue('L'), dims=2)
            # Change gap positions
            aln = shuffle_msa(msa, dims=2, fixedgaps=false)
            @test aln != getresidues(msa)
            @test (aln .== GAP) != gaps[i]
            @test lcol[i] == mean(aln .== Residue('L'), dims=1)
        end
    end

    @testset "Sequences" begin

        for  i in 1:N
            # Fixed gaps
            msa = msas[i]
            aln = shuffle_msa(msa, dims=1, fixedgaps=true)
            @test aln != getresidues(msa)
            @test (aln .== GAP) == gaps[i]
            @test lcol[i] != mean(aln .== Residue('L'), dims=1)
            @test lseq[i] == mean(aln .== Residue('L'), dims=2)
            # Change gap positions
            aln = shuffle_msa(msa, dims=1, fixedgaps=false)
            @test aln != getresidues(msa)
            @test (aln .== GAP) != gaps[i]
            @test lseq[i] == mean(aln .== Residue('L'), dims=2)
        end
    end

    @testset "Reference" begin
        
        for msa in msas
            ref = getsequence(msa, 1)
            for dims in [1, 2]
                for gaps in [true, false]
                    aln = shuffle_msa(msa, dims=dims, fixedgaps=gaps, fixed_reference=true)
                    @test getsequence(aln, 1) == ref
                end
            end
        end

        msa = msas[end] # AnnotatedMultipleSequenceAlignment
        aln = shuffle_msa(msa, dims=1, 1, fixed_reference=true)
        @test msa == aln
    end

    @testset "Subset" begin
        for msa in msas
            seqs_to_move = [3, 4]
            shuffled_seqs = shuffle_msa(msa, seqs_to_move, dims=1)
            @test msa[1:2, :] == shuffled_seqs[1:2, :]
            @test msa[seqs_to_move, :] != shuffled_seqs[seqs_to_move, :]

            cols_to_move = [9, 10, 11, 12]
            shuffled_cols = shuffle_msa(MersenneTwister(0), msa, cols_to_move, dims=2)
            @test msa[:, 1:8] == shuffled_cols[:, 1:8]
            @test msa[:, cols_to_move] != shuffled_cols[:, cols_to_move]

            # Annotations
            if isa(msa, AnnotatedMultipleSequenceAlignment)
                @test isempty(getannotsequence(shuffled_seqs, "C3N734_SULIY/1-95", "Shuffled", ""))
                @test isempty(getannotsequence(shuffled_seqs, "H2C869_9CREN/7-104", "Shuffled", ""))
                @test getannotsequence(shuffled_seqs, "Y070_ATV/2-70", "Shuffled", "") == "true"
                @test getannotsequence(shuffled_seqs, "F112_SSV1/3-112", "Shuffled", "") == "true"
                
                #                                                                         1111
                #                                                                1234567890123
                @test startswith(getannotcolumn(shuffled_cols, "Shuffled", ""), "0000000011110")

                # MIToS modifications
                any(startswith("2 sequences shuffled."), values(getannotfile(shuffled_seqs)))
                any(startswith("4 columns shuffled."), values(getannotfile(shuffled_cols)))
            end
        end

        # Shuffling a single sequence or column
        annot_msa = msas[end] # AnnotatedMultipleSequenceAlignment
        
        ref = getsequence(annot_msa, 1)
        @test getsequence(shuffle_msa(annot_msa, 1, dims=1), 1) != ref
        @test getsequence(shuffle_msa(annot_msa, "C3N734_SULIY/1-95", dims=1), 1) != ref

        # 10 == "15" : "EKQE"
        shuffled_col_a = shuffle_msa(annot_msa, 10, dims=2)[:, 10]
        shuffled_col_b = shuffle_msa(annot_msa, "15", dims=2)[:, "15"]
        @test replace(replace(join(shuffled_col_a), "E"=>""), "K" => "") == "Q"
        @test replace(replace(join(shuffled_col_b), "E"=>""), "K" => "") == "Q"
        @test annot_msa[:, 10] != shuffled_col_a
        @test annot_msa[:, "15"] != shuffled_col_b
    end
end
