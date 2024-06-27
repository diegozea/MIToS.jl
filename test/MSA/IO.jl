@testset "IO" begin

    msa_types = (
        Matrix{Residue},
        NamedResidueMatrix{Array{Residue,2}},
        MultipleSequenceAlignment,
        AnnotatedMultipleSequenceAlignment,
    )

    # Sequence from PF09645
    F112_SSV1 = collect(
        ".....QTLNSYKMAEIMYKILEKKGELTLEDILAQFEISVPSAYNIQRALKAICERHPD" *
        "ECEVQYKNRKTTFKWIKQEQKEEQKQEQTQDNIAKIFDAQPANFEQTDQGFIKAKQ.....",
    )

    # Test pfam stockholm parser using the 4 sequence full MSA for PF09645
    # > Order: Tree
    # > Inserts lower case
    # > Gaps as "." or "-" (mixed)
    # > Pfam version 28.0, based on UniProt release 2014_07
    @testset "Stockholm" begin

        pf09645_sto = joinpath(DATA, "PF09645_full.stockholm")
        rna_sto = joinpath(DATA, "upsk_rna.sto")
        clustal_sto = joinpath(DATA, "clustalo-I20240512-trunc.aln-stockholm")

        @testset "Read" begin

            # TODO : @inferred read_file(pf09645_sto, Stockholm);

            @test isa(read_file(pf09645_sto, Stockholm), AnnotatedMultipleSequenceAlignment)
            for T in msa_types
                @test isa(read_file(pf09645_sto, Stockholm, T), T)
            end

            for T in msa_types
                @test read_file(pf09645_sto, Stockholm, T) ==
                      read_file(pf09645_sto, Stockholm)
            end

            @testset "Empty lines" begin
                # This Rfam alignment has an empty line between the file annotations and the sequences
                for T in msa_types # NOTE: nucleotides as Residue objects
                    @test isa(read_file(rna_sto, Stockholm, T), T)
                end
            end
        end

        @testset "Output types" begin

            msa_objects = [read_file(pf09645_sto, Stockholm, T) for T in msa_types]

            @testset "Sequence Names" begin

                default = ["1", "2", "3", "4"]
                pfnames = [
                    "C3N734_SULIY/1-95",
                    "H2C869_9CREN/7-104",
                    "Y070_ATV/2-70",
                    "F112_SSV1/3-112",
                ]

                @test sequencenames(msa_objects[1]) == default
                for i = 2:4
                    @test sequencenames(msa_objects[i]) == pfnames
                end
            end

            @testset "Size" begin

                for msa in msa_objects
                    @test size(msa, 1) == 4
                    @test size(msa, 2) == sum(F112_SSV1 .!= Ref('.')) # without inserts
                    @test view(msa, 4, :) == map(Residue, F112_SSV1[F112_SSV1.!=Ref('.')])
                end

                msacl = read_file(
                    clustal_sto,
                    Stockholm,
                    generatemapping = true,
                    useidcoordinates = true,
                )
                @test size(msacl, 1) == 2
                @test maximum(getsequencemapping(msacl, "Q8BX79|reviewed|Probable")) == 349
                @test maximum(getsequencemapping(msacl, "Q923Y8|reviewed|Trace")) == 332
            end
        end

        @testset "Annotations" begin

            msa = read_file(pf09645_sto, Stockholm)

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

            @test parse_file(pfam_string, Stockholm) == read_file(pf09645_sto, Stockholm)

            msa = parse_file(pfam_string, Stockholm)
            @test !isempty(msa)
            @test length(annotations(msa)) > 8 # 8 + modifications
            @test length(getannotcolumn(msa)) == 2
            @test length(getannotresidue(msa)) == 1
        end

        @testset "File write Matrix{Residue}" begin

            path = tempdir()
            tmp_file = joinpath(path, ".tmp.stockholm")
            try
                msa = read_file(pf09645_sto, Stockholm, Matrix{Residue})
                write_file(tmp_file, msa, Stockholm)
                @test read_file(tmp_file, Stockholm, Matrix{Residue}) == msa
            finally
                if isfile(tmp_file)
                    rm(tmp_file)
                end
            end
        end

        @testset "Keep insert columns" begin

            msa = read_file(pf09645_sto, Stockholm, keepinserts = true)
            # Aligned columns
            @test (collect(getannotcolumn(msa, "Aligned")) .== Ref('1')) ==
                  (F112_SSV1 .!= Ref('.'))
            seq = "...mp---NSYQMAEIMYKILQQKKEISLEDILAQFEISASTAYNVQRTLRMICEKHPDECEVQTKNRRTIFKWIKNEETTEEGQEE--QEIEKILNAQPAE-------------k...."
            @test stringsequence(msa, 1) == replace(uppercase(seq), '.' => '-')

            @testset "Print inserts" begin

                io = IOBuffer()
                print_file(io, msa, Stockholm)
                printed = String(take!(io))
                @test occursin(seq, printed)
            end
        end
    end

    @testset "FASTA" begin

        pf09645_fas = joinpath(DATA, "PF09645_full.fasta.gz")
        gaoetal2011 = joinpath(DATA, "Gaoetal2011.fasta")

        @testset "Read" begin

            @test isa(read_file(pf09645_fas, FASTA), AnnotatedMultipleSequenceAlignment)
            for T in msa_types
                @test isa(read_file(pf09645_fas, FASTA, T), T)
            end

            for T in msa_types
                @test read_file(pf09645_fas, FASTA, T) == read_file(pf09645_fas, FASTA)
            end

            @testset "Download" begin

                @test read_file(gaoetal2011, FASTA) == read_file(
                    "https://raw.githubusercontent.com/diegozea/MIToS.jl/master/test/data/Gaoetal2011.fasta",
                    FASTA,
                )
                @test read_file(pf09645_fas, FASTA) == read_file(
                    "https://raw.githubusercontent.com/diegozea/MIToS.jl/master/test/data/PF09645_full.fasta.gz",
                    FASTA,
                )
            end
        end

        @testset "Output types" begin

            gaoetal_msas = [read_file(gaoetal2011, FASTA, T) for T in msa_types]
            pfam_msas = [read_file(pf09645_fas, FASTA, T) for T in msa_types]

            @testset "Sequence Names" begin

                @test sequencenames(gaoetal_msas[1]) == ["1", "2", "3", "4", "5", "6"]
                @test sequencenames(pfam_msas[1]) == ["1", "2", "3", "4"]
                for i = 2:4
                    @test sequencenames(gaoetal_msas[i]) ==
                          ["SEQ1", "SEQ2", "SEQ3", "SEQ4", "SEQ5", "SEQ6"]
                    @test sequencenames(pfam_msas[i]) == [
                        "C3N734_SULIY/1-95",
                        "H2C869_9CREN/7-104",
                        "Y070_ATV/2-70",
                        "F112_SSV1/3-112",
                    ]
                end
            end

            @testset "Residues" begin

                list_seq = [
                    res"DAWAEE",
                    res"DAWAEF",
                    res"DAWAED",
                    res"DAYCMD",
                    res"DAYCMT",
                    res"DAYCMT",
                ]
                residues = permutedims(hcat(list_seq...), [2, 1])

                for msa in gaoetal_msas

                    @test msa == residues
                    @test isa(getresidues(msa), Matrix{Residue})
                    @test getresidues(msa) == residues
                    @test getresiduesequences(msa) == list_seq
                    @test stringsequence(msa, 1) == "DAWAEE"
                end

                for msa in gaoetal_msas[2:end]
                    # getsequence returns a matrix, use vec(seq) or dropdims(seq,  dims=1)
                    # to get a vector:
                    @test vec(getresidues(getsequence(msa, 1))) == res"DAWAEE"
                    @test dropdims(getresidues(getsequence(msa, 1)), dims = 1) ==
                          res"DAWAEE"
                end
            end

            @testset "Size" begin

                for msa in gaoetal_msas
                    @test size(msa, 1) == 6
                    @test size(msa, 2) == 6
                    @test view(msa, 4, :) == res"DAYCMD"
                end

                for msa in pfam_msas
                    @test size(msa, 1) == 4
                    @test size(msa, 2) == sum(F112_SSV1 .!= Ref('.')) # without inserts
                    @test view(msa, 4, :) == map(Residue, F112_SSV1[F112_SSV1.!=Ref('.')])
                end
            end
        end

        @testset "String input/output" begin

            msa = read_file(gaoetal2011, FASTA)
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

            @test parse_file(fasta_string, FASTA) == msa

            out = IOBuffer()
            print_file(out, msa, FASTA)
            @test String(take!(out)) == fasta_string
        end

        @testset "File input/output" begin

            msa = read_file(gaoetal2011, FASTA)

            path = tempdir()
            uncompressed = joinpath(path, ".tmp.fasta")
            try
                write_file(uncompressed, msa, FASTA)
                @test read_file(uncompressed, FASTA) == msa
            finally
                if isfile(uncompressed)
                    rm(uncompressed)
                end
            end
            compressed = joinpath(path, ".tmp.fasta.gz")
            try
                write_file(compressed, msa, FASTA)
                @test read_file(compressed, FASTA) == msa
            finally
                if isfile(compressed)
                    rm(compressed)
                end
            end
        end

        @testset "Non standard residues and mapping" begin

            seqs =
                read_file(joinpath(DATA, "alphabet.fasta"), FASTA, generatemapping = true)

            @test vec(seqs[1, :]) == res"ARNDCQEGHILKMFPSTWYV"
            for i = 2:nsequences(seqs)
                @test vec(seqs[i, :]) == res"ARXDCQEGHILKMFPSTWYV"
                @test getsequencemapping(seqs, 1) == getsequencemapping(seqs, i)
            end
        end
    end

    @testset "Raw" begin

        # AnnotatedMultipleSequenceAlignment
        raw = read_file(joinpath(DATA, "gaps.txt"), Raw)
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
            @test stringsequence(raw, 10) == "T---------"

            for i = 1:10
                @test stringsequence(mat, i) == stringsequence(raw, i)
            end

            @test parse_file(raw_string, Raw) == raw
        end

        @testset "Print" begin

            buffer = IOBuffer()
            print_file(buffer, raw, Raw)
            @test String(take!(buffer)) == raw_string

            print_file(buffer, mat, Raw)
            @test String(take!(buffer)) == raw_string
        end

        @testset "Stats" begin

            @test gapfraction(mat) == 0.45
            @test vec(gapfraction(mat, 1)) ≈
                  [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
            @test vec(gapfraction(mat, 2)) ≈
                  [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

            @test residuefraction(mat) == 0.55
            @test vec(residuefraction(mat, 1)) ≈
                  [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
            @test vec(residuefraction(mat, 2)) ≈
                  [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]

            @test vec(coverage(mat)) ≈ [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]

            @test vec(columngapfraction(mat)) ≈
                  [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
        end

        @testset "Reference and gapstrip" begin

            gs = gapstrip(mat, coveragelimit = 0.5, gaplimit = 0.5)
            @test vec(getsequence(gs, 1)) == res"THAYQAIH"
            @test ncolumns(gs) == 8
            @test nsequences(gs) == 6

            ref = setreference!(copy(mat), 2)
            @test vec(getsequence(ref, 1)) == res"THAYQAIHQ-"
            ref = adjustreference(ref)
            @test vec(getsequence(ref, 1)) == res"THAYQAIHQ"
        end
    end

    @testset "NBRF/PIR" begin

        # Example NBRF file: http://iubio.bio.indiana.edu/soft/molbio/readseq/classic/src/Formats
        # "The sequence is free format and may be interrupted by blanks for ease of reading"
        example = joinpath(DATA, "example.nbrf")
        # Alignment file (PIR) from https://salilab.org/modeller/9v7/manual/node445.html
        modeller = joinpath(DATA, "modeller.pir.gz")
        # http://emboss.sourceforge.net/docs/themes/seqformats/NbrfFormat.html
        # "sequence may contain punctuation symbols to indicate various degrees of
        # reliability of the data"
        emboss = joinpath(DATA, "emboss.pir")
        # http://caps.ncbs.res.in/pass2v3/pir.html
        # Example from pass2 with spaces added at the end of the id lines.
        pass2 = joinpath(DATA, "pass2.pir")

        @testset "Read" begin

            @test isa(read_file(modeller, PIR), AnnotatedMultipleSequenceAlignment)
            for T in msa_types
                @test isa(read_file(modeller, PIR, T), T)
            end

            msa = read_file(modeller, PIR)
            @test size(msa) == (2, 106)
        end

        @testset "Print" begin

            msa = read_file(example, PIR)

            @test stringsequence(msa, 1) ==
                  "MTNIRKSHPLFKIINHSFIDLPAPSVTHICRDVNYGWLIRYTWIGGQPVEHPFIIIGQLASISYFSIILILMPISGIVEDKMLKWN"

            out = IOBuffer()
            print_file(out, msa, PIR)
            @test String(take!(out)) == """
                >P1;CBRT
                Cytochrome b - Rat mitochondrion (SGC1)
                MTNIRKSHPLFKIINHSFIDLPAPSVTHICRDVNYGWLIRYTWIGGQPVEHPFIIIGQLASISYFSIILILMPISGIVED
                KMLKWN*
                """

            out = IOBuffer()
            print_file(out, msa.matrix, PIR) # NamedResidueMatrix{Array{Residue,2}}
            @test String(take!(out)) == """
                >XX;CBRT

                MTNIRKSHPLFKIINHSFIDLPAPSVTHICRDVNYGWLIRYTWIGGQPVEHPFIIIGQLASISYFSIILILMPISGIVED
                KMLKWN*
                """

            out = IOBuffer()
            print_file(out, msa.matrix.array, PIR) # Matrix{Residue}
            @test String(take!(out)) == """
                >XX;1

                MTNIRKSHPLFKIINHSFIDLPAPSVTHICRDVNYGWLIRYTWIGGQPVEHPFIIIGQLASISYFSIILILMPISGIVED
                KMLKWN*
                """

        end

        @testset "Gaps" begin

            msa = read_file(modeller, PIR)
            @test gapfraction(msa, 2) ≈ [0.0, 0.49056603773584906]
            @test getannotsequence(msa, "5fd1", "Type") == "P1"
            @test getannotsequence(msa, "5fd1", "Title") ==
                  "structureX:5fd1:1    :A:106  :A:ferredoxin:Azotobacter vinelandii: 1.90: 0.19"
        end

        @testset "String input/output" begin

            emboss_string = """
                >P1;AZBR
                finger protein zfpA - turnip fern chloroplast
                GDVE(G.K.G.I.F=T,M,C.S.Q,C.H.V,E.K.G.G.K.H)
                FTGPNLHGLFGRK.TGQAVGYSYTAANK.NK.GIIWGDDTLM
                EYLENPK.RYIPGTK.MVFTGLSK.YRE
                RTNLIAYLK.EK.TAA*
                """
            # MIToS always reads . as -
            msa = read_file(emboss, PIR, deletefullgaps = false)
            @test parse_file(emboss_string, PIR, deletefullgaps = false) == msa

            out = IOBuffer()
            print_file(out, msa, PIR)
            @test String(take!(out)) == """
                >P1;AZBR
                finger protein zfpA - turnip fern chloroplast
                GDVEG-K-G-I-FTMC-S-QC-H-VE-K-G-G-K-HFTGPNLHGLFGRK-TGQAVGYSYTAANK-NK-GIIWGDDTLMEY
                LENPK-RYIPGTK-MVFTGLSK-YRERTNLIAYLK-EK-TAA*
                """
        end

        @testset "Spaces at the end of the id line" begin
            lines = readlines(pass2)
            @test lines[1] == ">P1;1bbha- "
            @test lines[5] == ">P1;1cpq-- "
            @test lines[9] == ">P1;256bb- "
            msa = read_file(pass2, PIR)
            @test sequencenames(msa) == ["1bbha-", "1cpq--", "256bb-"]
        end

        @testset "Duplicated identifiers" begin

            duplicated_ids_file = joinpath(DATA, "duplicated_ids.pir")
            @test_throws ArgumentError read_file(duplicated_ids_file, PIR)
        end
    end

    @testset "A3M" begin

        @testset "Adding insert gaps" begin
            # Simple tests for adding insert gaps in different locations
            @test MSA._add_insert_gaps!(["MAHI", "MgAHI"]) == ["M.AHI", "MgAHI"]
            @test MSA._add_insert_gaps!(["MAHILLI", "MAHIgLLIhhh", "MAHILLI", "MAHILLI"]) ==
                  ["MAHI.LLI...", "MAHIgLLIhhh", "MAHI.LLI...", "MAHI.LLI..."]
            @test MSA._add_insert_gaps!(["MAhI", "MALh"]) == ["MAhI.", "MA.Lh"]
            @test MSA._add_insert_gaps!(["MALI", "hhhMALI"]) == ["...MALI", "hhhMALI"]
            @test MSA._add_insert_gaps!(["MAggLLI", "MAgggLLI", "MAgLLI", "MALLI"]) ==
                  ["MAgg.LLI", "MAgggLLI", "MAg..LLI", "MA...LLI"]
        end

        @testset "Single insert column" begin
            # A3M example from https://yanglab.qd.sdu.edu.cn/trRosetta/msa_format.html
            #                                 |insert
            seq = "-----RTKRLREAVRVYLAENGrSHTVDIFDHLNDRFSWGATMNQVGNILAKDNRFEKVGHVRD-FFRGARYTVCVWDLAS-----------"
            seq_without_insert = replace(seq, "r" => "")

            msa = read_file(joinpath(DATA, "yanglab.a3m"), A3M)

            @test ncolumns(msa) == 91 # 92 with the insert column
            @test nsequences(msa) == 7
            @test stringsequence(msa, "6") == seq_without_insert

            @testset "Print inserts" begin

                @testset "There are no inserts" begin
                    io = IOBuffer()
                    print_file(io, msa, A3M)
                    printed = String(take!(io))
                    # The insert column has been removed when reading the MSA
                    @test occursin(seq_without_insert, printed)
                end

                @testset "Keep inserts" begin
                    # Keeping the insert column when reading the MSA
                    msa = read_file(joinpath(DATA, "yanglab.a3m"), A3M, keepinserts = true)
                    @test ncolumns(msa) == 92

                    io = IOBuffer()
                    print_file(io, msa, A3M)
                    printed = String(take!(io))
                    @test occursin(seq, printed)
                    # Test that gaps are not added to the insert column for A3M
                    @test !occursin(".", printed)

                    @testset "Saving as A2M" begin
                        io = IOBuffer()
                        print_file(io, msa, A2M)
                        printed = String(take!(io))
                        # A2M keeps the insert column
                        @test occursin(seq, printed)
                        # and add explicit gaps for inserts
                        @test occursin(".", printed)
                        # Ensure the correct location of the gaps
                        @test occursin("RP.RN", printed)
                    end
                end
            end
        end
    end
end
