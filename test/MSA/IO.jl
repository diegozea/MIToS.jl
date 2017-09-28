@testset "IO" begin

    msa_types = (
        Matrix{Residue},
        NamedArray{Residue,2},
        MultipleSequenceAlignment,
        AnnotatedMultipleSequenceAlignment
        )

    # Sequence from PF09645
    F112_SSV1 = collect(".....QTLNSYKMAEIMYKILEKKGELTLEDILAQFEISVPSAYNIQRALKAICERHPD" *
                        "ECEVQYKNRKTTFKWIKQEQKEEQKQEQTQDNIAKIFDAQPANFEQTDQGFIKAKQ.....")

    # Test pfam stockholm parser using the 4 sequence full MSA for PF09645
    # > Order: Tree
    # > Inserts lower case
    # > Gaps as "." or "-" (mixed)
    # > Pfam version 28.0, based on UniProt release 2014_07
    @testset "Stockholm" begin

        pf09645_sto = joinpath(pwd(), "data", "PF09645_full.stockholm")

        @testset "Read" begin

            @test isa(read(pf09645_sto, Stockholm), AnnotatedMultipleSequenceAlignment)
            for T in msa_types
                @test isa(read(pf09645_sto, Stockholm, T), T)
            end

            for T in msa_types
                @test read(pf09645_sto, Stockholm, T) == read(pf09645_sto, Stockholm)
            end
        end

        @testset "Output types" begin

            msa_objects = [ read(pf09645_sto, Stockholm, T)  for T in msa_types ]

            @testset "Sequence Names" begin

                default = ["1", "2", "3", "4"]
                pfnames = [ "C3N734_SULIY/1-95","H2C869_9CREN/7-104",
                            "Y070_ATV/2-70",    "F112_SSV1/3-112"   ]

                @test sequencenames(msa_objects[1]) == default
                for i in 2:4
                    @test sequencenames(msa_objects[i]) == pfnames
                end
            end

            @testset "Size" begin

                for msa in msa_objects
                    @test size(msa, 1) == 4
                    @test size(msa, 2) == sum(F112_SSV1 .!= '.') # without inserts
                    @test view(msa,4,:) == map(Residue, F112_SSV1[ F112_SSV1 .!= '.' ])
                end
            end
        end

        @testset "Annotations" begin

            msa = read(pf09645_sto, Stockholm)

            @test !isempty(msa)
            @test length(annotations(msa)) > 8 # 8 + modifications
            @test length(getannotcolumn(msa)) == 2
            @test length(getannotresidue(msa)) == 1
            @test getannotresidue(msa, "F112_SSV1/3-112", "SS") ==
                "X---HHHHHHHHHHHHHHHSEE-HHHHHHHH---HHHHHHHHHHHHHHHHH-TTTEEEEE-SS-EEEEE--XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
            @test getannotcolumn(msa, "seq_cons") ==
                "...NshphAclhaKILppKtElolEDIlAQFEISsosAYsI.+sL+hICEpH.-ECpsppKsRKTlhh.hKpEphppptpEp..ppItKIhsAp................"
            @test getannotsequence(msa, "F112_SSV1/3-112", "DR") == "PDB; 2VQC A; 4-73;"
        end

        @testset "String input" begin

            pfam_string = """
                # STOCKHOLM 1.0
                #=GS C3N734_SULIY/1-95   AC C3N734.1
                #=GS H2C869_9CREN/7-104  AC H2C869.1
                #=GS Y070_ATV/2-70       AC Q3V4T1.1
                #=GS F112_SSV1/3-112     AC P20220.1
                #=GS F112_SSV1/3-112     DR PDB; 2VQC A; 4-73;
                C3N734_SULIY/1-95              ...mp---NSYQMAEIMYKILQQKKEISLEDILAQFEISASTAYNVQRTLRMICEKHPDECEVQTKNRRTIFKWIKNEETTEEGQEE--QEIEKILNAQPAE-------------k....
                H2C869_9CREN/7-104             ...nk--LNDVQRAKLLVKILQAKGELDVYDIMLQFEISYTRAIPIMKLTRKICEAQ-EICTYDEKEHKLVSLNAKKEKVEQDEEENEREEIEKILDAH----------------trreq
                Y070_ATV/2-70                  qsvne-------VAQQLFSKLREKKEITAEDIIAIYNVTPSVAYAIFTVLKVMCQQHQGECQAIKRGRKTVI-------------------------------------------vskq.
                F112_SSV1/3-112                .....QTLNSYKMAEIMYKILEKKGELTLEDILAQFEISVPSAYNIQRALKAICERHPDECEVQYKNRKTTFKWIKQEQKEEQKQEQTQDNIAKIFDAQPANFEQTDQGFIKAKQ.....
                #=GR F112_SSV1/3-112     SS    .....X---HHHHHHHHHHHHHHHSEE-HHHHHHHH---HHHHHHHHHHHHHHHHH-TTTEEEEE-SS-EEEEE--XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX.....
                #=GC SS_cons                   .....X---HHHHHHHHHHHHHHHSEE-HHHHHHHH---HHHHHHHHHHHHHHHHH-TTTEEEEE-SS-EEEEE--XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX.....
                #=GC seq_cons                  ........NshphAclhaKILppKtElolEDIlAQFEISsosAYsI.+sL+hICEpH.-ECpsppKsRKTlhh.hKpEphppptpEp..ppItKIhsAp................h....
                //
                """

            @test parse(pfam_string, Stockholm) == read(pf09645_sto, Stockholm)

            msa = parse(pfam_string, Stockholm)
            @test !isempty(msa)
            @test length(annotations(msa)) > 8 # 8 + modifications
            @test length(getannotcolumn(msa)) == 2
            @test length(getannotresidue(msa)) == 1
        end

        @testset "Keep insert columns" begin

            msa = read(pf09645_sto, Stockholm, keepinserts=true)
            # Aligned columns
            @test (collect(getannotcolumn(msa, "Aligned")) .== '1') == (F112_SSV1 .!= '.')
            @test stringsequence(msa, 1) == replace(uppercase(
                "...mp---NSYQMAEIMYKILQQKKEISLEDILAQFEISASTAYNVQRTLRMICEKHPDECEVQTKNRRTIFKWIKNEETTEEGQEE--QEIEKILNAQPAE-------------k...."
                ), '.', '-')
        end
    end

    @testset "FASTA" begin

        pf09645_fas = joinpath(pwd(), "data", "PF09645_full.fasta.gz")
        gaoetal2011 = joinpath(pwd(), "data", "Gaoetal2011.fasta")

        @testset "Read" begin

            @test isa(read(pf09645_fas, FASTA), AnnotatedMultipleSequenceAlignment)
            for T in msa_types
                @test isa(read(pf09645_fas, FASTA, T), T)
            end

            for T in msa_types
                @test read(pf09645_fas, FASTA, T) == read(pf09645_fas, FASTA)
            end

            @testset "Download" begin

                @test read(gaoetal2011, FASTA) ==
                    read("https://raw.githubusercontent.com/diegozea/MIToS.jl/master/test/data/Gaoetal2011.fasta",
                        FASTA)
                @test read(pf09645_fas, FASTA) ==
                    read("https://raw.githubusercontent.com/diegozea/MIToS.jl/master/test/data/PF09645_full.fasta.gz",
                        FASTA)
            end
        end

        @testset "Output types" begin

            gaoetal_msas = [ read(gaoetal2011, FASTA, T)  for T in msa_types ]
            pfam_msas  = [ read(pf09645_fas, FASTA, T)  for T in msa_types ]

            @testset "Sequence Names" begin

                @test sequencenames(gaoetal_msas[1]) == ["1", "2", "3", "4", "5", "6"]
                @test sequencenames(pfam_msas[1])  == ["1", "2", "3", "4"]
                for i in 2:4
                    @test sequencenames(gaoetal_msas[i]) ==
                        ["SEQ1", "SEQ2", "SEQ3", "SEQ4", "SEQ5", "SEQ6"]
                    @test sequencenames(pfam_msas[i]) == [  "C3N734_SULIY/1-95",
                                                            "H2C869_9CREN/7-104",
                                                            "Y070_ATV/2-70",
                                                            "F112_SSV1/3-112"   ]
                end
            end

            @testset "Residues" begin

                list_seq = [res"DAWAEE",
                            res"DAWAEF",
                            res"DAWAED",
                            res"DAYCMD",
                            res"DAYCMT",
                            res"DAYCMT" ]
                residues = permutedims(hcat(list_seq...), [2,1])

                for msa in gaoetal_msas

                    @test msa == residues
                    @test isa(getresidues(msa), Matrix{Residue})
                    @test getresidues(msa) == residues
                    @test getresiduesequences(msa) == list_seq
                    @test stringsequence(msa, 1) == "DAWAEE"
                end

                for msa in gaoetal_msas[2:end]
                    # getsequence returns a matrix, use vec(seq) or squeeze(seq,1)
                    # to get a vector:
                    @test vec(getresidues(getsequence(msa, 1))) == res"DAWAEE"
                    @test squeeze(getresidues(getsequence(msa, 1)),1) == res"DAWAEE"
                end
            end

            @testset "Size" begin

                for msa in gaoetal_msas
                    @test size(msa, 1) == 6
                    @test size(msa, 2) == 6
                    @test view(msa,4,:) == res"DAYCMD"
                end

                for msa in pfam_msas
                    @test size(msa, 1) == 4
                    @test size(msa, 2) == sum(F112_SSV1 .!= '.') # without inserts
                    @test view(msa,4,:) == map(Residue, F112_SSV1[ F112_SSV1 .!= '.' ])
                end
            end
        end

        @testset "String input/output" begin

            msa = read(gaoetal2011, FASTA)
            fasta_string = """
                >SEQ1
                DAWAEE
                >SEQ2
                DAWAEF
                >SEQ3
                DAWAED
                >SEQ4
                DAYCMD
                >SEQ5
                DAYCMT
                >SEQ6
                DAYCMT
                """

            @test parse(fasta_string, FASTA) == msa

            out = IOBuffer()
            print(out, msa, FASTA)
            @test takebuf_string(out) == fasta_string
        end

        @testset "File input/output" begin

            msa = read(gaoetal2011, FASTA)

            path = tempdir()
            uncompressed = joinpath(path, ".tmp.fasta")
            try
                write(uncompressed, msa, FASTA)
                @test read(uncompressed, FASTA) == msa
            finally
                if isfile(uncompressed)
                    rm(uncompressed)
                end
            end
            compressed = joinpath(path, ".tmp.fasta.gz")
            try
                write(compressed, msa, FASTA)
                @test read(compressed, FASTA) == msa
            finally
                if isfile(compressed)
                    rm(compressed)
                end
            end
        end

        @testset "Non standard residues and mapping" begin

            seqs = read(joinpath(pwd(), "data", "alphabet.fasta"), FASTA,
                        generatemapping=true)

            @test vec(seqs[1,:]) == res"ARNDCQEGHILKMFPSTWYV"
            for i in 2:nsequences(seqs)
                @test vec(seqs[i,:]) == res"ARXDCQEGHILKMFPSTWYV"
                @test getsequencemapping(seqs,1) == getsequencemapping(seqs,i)
            end
        end
    end

    @testset "Raw" begin

        # AnnotatedMultipleSequenceAlignment
        raw = read(joinpath(pwd(), "data", "gaps.txt"), Raw)
        mat = getresidues(raw)
        raw_string = """THAYQAIHQV
                        THAYQAIHQ-
                        THAYQAIH--
                        THAYQAI---
                        THAYQA----
                        THAYQ-----
                        THAY------
                        THA-------
                        TH--------
                        T---------
                        """
        @testset "Parse" begin

            @test stringsequence(raw, "1") == "THAYQAIHQV"
            @test stringsequence(raw, 10)  == "T---------"

            for i in 1:10
                @test stringsequence(mat, i)  == stringsequence(raw, i)
            end

            @test parse(raw_string, Raw) == raw
        end

        @testset "Print" begin

            buffer = IOBuffer()
            print(buffer, raw, Raw)
            @test takebuf_string(buffer) == raw_string

            print(buffer, mat, Raw)
            @test takebuf_string(buffer) == raw_string
        end

        @testset "Stats" begin

            @test gapfraction(mat) == 0.45
            @test vec(gapfraction(mat, 1)) ≈ [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
            @test vec(gapfraction(mat, 2)) ≈ [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]

            @test residuefraction(mat) == 0.55
            @test vec(residuefraction(mat, 1)) ≈ [1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1]
            @test vec(residuefraction(mat, 2)) ≈ [1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1]

            @test vec(coverage(mat)) ≈ [1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1]

            @test vec(columngapfraction(mat)) ≈ [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
        end

        @testset "Reference and gapstrip" begin

            gs = gapstrip(mat, coveragelimit=0.5, gaplimit=0.5)
            @test vec(getsequence(gs, 1)) == res"THAYQAIH"
            @test ncolumns(gs) == 8
            @test nsequences(gs) == 6

            ref = setreference!(copy(mat), 2)
            @test vec(getsequence(ref, 1)) == res"THAYQAIHQ-"
            ref = adjustreference(ref)
            @test vec(getsequence(ref, 1)) == res"THAYQAIHQ"
        end
    end
end
