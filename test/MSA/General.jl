@testset "Test using MSA files" begin

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

    @testset "getindex" begin

        residues = permutedims(hcat(res"DAWAEE",
                                    res"DAWAEF",
                                    res"DAWAED",
                                    res"DAYCMD",
                                    res"DAYCMT",
                                    res"DAYCMT"), [2,1])

        for index in [ 2, (2,:), (:,2), (2,2), (2:4, 2:4), (2,[3,4,5]),
                      (:, [1,2,3,4,5,6] .< 4 ), ([1,2,3,4,5,6] .< 4, :),
                      (:,:), ([1,2,3,4,5,6] .< 4, [1,2,3,4,5,6] .< 4) ]

            for msa in gaoetal_msas
                @test msa[index...] == residues[index...]
            end
        end
    end

    @testset "setindex! and copy" begin

        for msa in gaoetal_msas
            copy_msa = copy(msa)
            deepcopy_msa = deepcopy(msa)

            for (index, value) in [ ((1), Residue('H')),
                                    ((:,1), res"HHHHHH") ]

                deepcopy_msa[index...] = value
                @test deepcopy_msa[index...] == value
                @test msa[index...] != value

                copy_msa[index...] = value
                @test copy_msa[index...] == value
                @test msa[index...] != value
            end

            for (index, value) in [ ( 1 , Residue('H')),
                                    ( 4 , Residue('X')),
                                    ( : , res"HHHHHH" ) ]

                seq = copy(getsequence(msa, 4))
                seq[1, index] = value
                @test seq[1, index] == value
                @test msa[4, index] != value # since seq is a copy
            end
        end
    end

    @testset "Size" begin

        for aln in gaoetal_msas
            @test size(aln) == (6,6)
            @test length(aln) == 36
            @test ncolumns(aln) == 6
            @test nsequences(aln) == 6
            @test ncolumns(getsequence(aln, 4)) == 6
        end

        for aln in pfam_msas
            @test size(aln) == (4,110)
            @test length(aln) == 440
            @test ncolumns(aln) == 110
            @test nsequences(aln) == 4
        end
    end

    @testset "AnnotatedAlignedSequence and AlignedSequence" begin

        @testset "Creation" begin
            seq_types = (   Matrix{Residue},
                            NamedResidueMatrix{Array{Residue,2}},
                            AlignedSequence,
                            AnnotatedAlignedSequence    )

            for i in 1:length(msa_types)
                M = msa_types[i]
                S = seq_types[i]
                msa = pfam_msas[i]

                if M != Matrix{Residue}
                    for id in [ "C3N734_SULIY/1-95","H2C869_9CREN/7-104",
                                "Y070_ATV/2-70",    "F112_SSV1/3-112"   ]
                        annseq = getsequence(msa, id)
                        @test msa[id,:] == vec(annseq) # Sequences are matrices
                        @test isa(annseq, S)
                    end
                end

                for seq in 1:4
                    annseq = getsequence(msa, seq)
                    @test msa[seq,:] == vec(annseq) # Sequences are matrices
                    @test isa(annseq, S)
                end
            end
        end

        @testset "Annotations" begin
            msa = read(pf09645_sto, Stockholm)

            @test   getannotcolumn(msa, "SS_cons") ==
                    getannotcolumn(getsequence(msa, 4), "SS_cons")
            @test   getannotresidue(msa, "F112_SSV1/3-112", "SS") ==
                    getannotresidue(getsequence(msa, 4), "F112_SSV1/3-112", "SS")
            @test   getannotfile(msa) == getannotfile(getsequence(msa, 4))
            @test   getannotsequence(msa,"F112_SSV1/3-112","DR") ==
                    getannotsequence(getsequence(msa, 4), "F112_SSV1/3-112", "DR")
        end
    end

    @testset "Print" begin
        pfam = read(pf09645_sto, Stockholm)
        gao  = read(gaoetal2011, FASTA)

        @test stringsequence(gao, 4) == "DAYCMD"
        @test stringsequence(pfam, 1) == stringsequence(pfam, "C3N734_SULIY/1-95")

        for T in (Stockholm, FASTA, Raw)
            buffer = IOBuffer()
            print(buffer, pfam, T)
            @test parse(String(take!(buffer)), T) == pfam

            print(buffer, gao, T)
            @test parse(String(take!(buffer)), T) == gao
        end
    end
end
