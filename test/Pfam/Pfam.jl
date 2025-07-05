@testset "Download" begin

    pfam_code = "PF11591"
    @test_throws ErrorException downloadpfam("2vqc")
    filename = downloadpfam(pfam_code, filename = tempname() * ".gz")
    try
        aln = read_file(filename, Stockholm)
        if size(aln) == (6, 34)
            @test getannotfile(aln, "ID") == "2Fe-2S_Ferredox"
        end
    finally
        rm(filename)
    end
end

@testset "PDB code" begin

    msa = read_file(joinpath(DATA, "PF09645_full.stockholm"), Stockholm)
    @test getseq2pdb(msa)["F112_SSV1/3-112"] == [("2VQC", "A")]
end

@testset "Mapping PDB/Pfam" begin

    msa_file = joinpath(DATA, "PF09645_full.stockholm")
    sifts_file = joinpath(DATA, "2vqc.xml")
    pdb_file = joinpath(DATA, "2VQC.xml")
    msa = read_file(msa_file, Stockholm, generatemapping = true, useidcoordinates = true)
    cmap = msacolumn2pdbresidue(msa, "F112_SSV1/3-112", "2VQC", "A", "PF09645", sifts_file)
    res = residuesdict(read_file(pdb_file, PDBML), model = "1", chain = "A", group = "ATOM")

    #     -45              20 pdb
    #.....QTLNSYKMAEIMYKILEK  msa seq
    #     123456789012345678  msa col
    #     345678901234567890  uniprot 3-20
    #    ****              *
    #12345678901234567890123  ColMap

    @test_throws KeyError cmap[5]  # insert
    @test cmap[6] == ""           # missing
    @test cmap[7] == "4"
    @test cmap[8] == "5"
    @test cmap[23] == "20"

    #.....QTLNSYKMAEIMYKILEKKGELTLEDILAQFEISVPSAYNIQRALKAICERHPDECEVQYKNRKTTFKWIKQEQKEEQKQEQTQDNIAKIFDAQPANFEQTDQGFIKAKQ..... msa seq
    #.....X---HHHHHHHHHHHHHHHSEE-HHHHHHHH---HHHHHHHHHHHHHHHHH-TTTEEEEE-SS-EEEEE--XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX..... pdb ss
    #                                                                                                   111111111111111111111 ColMap hundreds
    #         111111111122222222223333333333444444444455555555556666666666777777777788888888889999999999000000000011111111112 ColMap tens
    #123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890 ColMap ones
    #                                                                           **                                     **

    @test_throws KeyError cmap[116]   # insert
    @test cmap[115] == ""             # missing
    @test cmap[77] == ""             # missing
    @test cmap[76] == "73"

    @testset "Residues" begin

        msares = msaresidues(msa, res, cmap)

        @test_throws KeyError msares[5]   # insert
        @test_throws KeyError msares[6]   # missing
        @test msares[7].id.name == "THR" # T
        @test msares[8].id.name == "LEU" # L
        @test msares[23].id.name == "LYS" # K

        @test_throws KeyError msares[116]   # insert
        @test_throws KeyError msares[117]   # missing
        @test_throws KeyError msares[77]    # missing
        @test msares[76].id.name == "LYS"   # K
    end

    @testset "Contacts" begin

        @test msacolumn2pdbresidue(
            msa,
            "F112_SSV1/3-112",
            "2VQC",
            "A",
            "PF09645",
            sifts_file,
            strict = true,
            checkpdbname = true,
        ) == cmap

        @test length(res) == 70
        @test length(unique(values(cmap))) == 71

        #     -45              20 pdb
        #.....QTLNSYKMAEIMYKILEK  msa seq
        #     123456789012345678  msa col
        #12345678901234567890123  ColMap
        #     345678901234567890  uniprot 3-20
        #     ***              *

        @test_throws KeyError res["3"]
        @test res[cmap[7]].id.number == "4"

        contacts = msacontacts(msa, res, cmap, 6.05)
        missings = sum(Int.(isnan.(contacts)), dims = 1)

        @test size(contacts) == (110, 110)

        @test missings[1] == 110              # missing
        @test missings[2] == (110 - (70 - 1)) # 70 residues, 1 diagonal NaN
        @test missings[3] == (110 - (70 - 1))
        @test missings[18] == (110 - (70 - 1))

        @test missings[110] == 110                 # missing
        @test missings[72] == 110                 # missing
        @test missings[71] == (110 - (70 - 1))

        ncontacts = sum(Int.(contacts .== 1.0), dims = 1)

        @test ncontacts[1] == 0
        @test ncontacts[2] == 2
        @test ncontacts[3] == 6

        @testset "using msaresidues" begin
            # Test MSA contact map using PDBResidues from the MSA

            msares = msaresidues(msa, res, cmap)
            ncol = ncolumns(msa)
            colmap = getcolumnmapping(msa)

            for i = 1:(ncol-1), j = (i+1):ncol
                if (colmap[i] in keys(msares)) && (colmap[j] in keys(msares))
                    @test (contacts[i, j] .== 1.0) ==
                          contact(msares[colmap[i]], msares[colmap[j]], 6.05)
                end
            end

            @testset "hasresidues" begin

                mask = hasresidues(msa, cmap)

                @test mask[1] == false
                @test mask[2] == true
                @test sum(mask) == length(msares)
            end
        end

        @testset "AUC" begin
            auc_value = AUC(
                buslje09(msa, lambda = 0.05, threshold = 62.0, samples = 0)[2],
                contacts,
            )
            @test isapprox(auc_value, 0.5291; atol = 3e-4)
        end
    end
end

@testset "AUC" begin

    ntru = 90
    nfal = 100
    score_tru = Float16[2 + 2x for x in randn(ntru)]
    score_fal = Float16[-2 + 2x for x in randn(nfal)]
    msacontacts = NamedArray(
        PairwiseListMatrix{Float16,false,Vector{Float16}}(
            vcat(ones(Float16, ntru), zeros(Float16, nfal)),
            Float16[1.0 for x = 1:20],
            20,
        ),
    )
    score = NamedArray(
        PairwiseListMatrix{Float16,false,Vector{Float16}}(
            vcat(score_tru, score_fal),
            Float16[NaN for x = 1:20],
            20,
        ),
    )
    correct = 1.0 - auc(roc(score_tru, score_fal))

    @test AUC(score, msacontacts) == correct

    score[1:end, 2] .= NaN
    msacontacts[1:end, 3] .= NaN
    list_score = getlist(getarray(score))
    list_contact = getlist(getarray(msacontacts))
    n_values = length(list_score)
    @test n_values == length(list_contact)
    tar = [
        list_score[i] for i = 1:n_values if
        !isnan(list_score[i]) & !isnan(list_score[i]) & (list_contact[i] == 1.0)
    ]
    non = [
        list_score[i] for i = 1:n_values if
        !isnan(list_score[i]) & !isnan(list_score[i]) & (list_contact[i] == 0.0)
    ]
    @test AUC(score, msacontacts) == AUC(roc(tar, non))
end

@testset "PDB mapping branches" begin
    msa_base = read_file(
        joinpath(DATA, "PF09645_full.stockholm"),
        Stockholm,
        generatemapping = true,
        useidcoordinates = true,
    )


    @testset "getseq2pdb duplicate ID" begin
        msa_dup = read_file(joinpath(DATA, "PF09645_full.stockholm"), Stockholm)
        setannotsequence!(msa_dup, "F112_SSV1/3-112", "DR", "PDB; 9XYZ B; 10-20;")
        seq2pdb = getseq2pdb(msa_dup)
        @test seq2pdb["F112_SSV1/3-112"] == [("2VQC", "A"), ("9XYZ", "B")]
    end

    @testset "skip blank Pfnum" begin
        # residue 2 in dummy_sifts__file.xml has an empty Pfam number so it
        # should be ignored in the mapping
        xml_path = joinpath(DATA, "dummy_sifts__file.xml")
        map_skip = msacolumn2pdbresidue(
            msa_base,
            "F112_SSV1/3-112",
            "2VQC",
            "A",
            "PF09645",
            xml_path;
            strict = false,
            checkpdbname = false,
            missings = true,
        )
        @test length(map_skip) == 1
        @test minimum(keys(map_skip)) == 6
    end

    @testset "warning vs strict" begin
        # Pfam number 4 maps to PDB residue SER while the MSA has T.
        # With strict=false a warning is expected; strict=true throws an error.
        xml_path = joinpath(DATA, "dummy_sifts__file.xml")
        @test_logs (:warn, r"MSA sequence residue") match_mode = :any begin
            msacolumn2pdbresidue(
                msa_base,
                "F112_SSV1/3-112",
                "2VQC",
                "A",
                "PF09645",
                xml_path;
                strict = false,
                checkpdbname = false,
                missings = true,
            )
        end
        @test_throws ErrorException msacolumn2pdbresidue(
            msa_base,
            "F112_SSV1/3-112",
            "2VQC",
            "A",
            "PF09645",
            xml_path;
            strict = true,
            checkpdbname = false,
            missings = true,
        )
    end

    @testset "PDB name check" begin
        # The PDB residue at Pfam number 4 is SER but the MSA residue is T, so
        # enabling checkpdbname should raise an error regardless of strict mode.
        xml_path = joinpath(DATA, "dummy_sifts__file.xml")
        @test_throws ErrorException msacolumn2pdbresidue(
            msa_base,
            "F112_SSV1/3-112",
            "2VQC",
            "A",
            "PF09645",
            xml_path;
            strict = false,
            checkpdbname = true,
            missings = true,
        )
    end

    @testset "overload equivalence" begin
        
        mktempdir() do tmpfolder
            cd(tmpfolder) do 
                # Compare the results from the different overloads.  Two of them
                # download the SIFTS file on demand.
                setannotfile!(msa_base, "AC", "PF09645.1")
                result1 = msacolumn2pdbresidue(
                    msa_base,
                    "F112_SSV1/3-112",
                    "2VQC",
                    "A",
                    "PF09645";
                    strict = false,
                    checkpdbname = false,
                    missings = true,
                )
                result2 = msacolumn2pdbresidue(
                    msa_base,
                    "F112_SSV1/3-112",
                    "2VQC",
                    "A",
                    "PF09645",
                    joinpath(DATA, "2vqc.xml");
                    strict = false,
                    checkpdbname = false,
                    missings = true,
                )
                result3 = msacolumn2pdbresidue(
                    msa_base,
                    "F112_SSV1/3-112",
                    "2VQC",
                    "A";
                    strict = false,
                    checkpdbname = false,
                    missings = true,
                )
                @test result1 == result2 == result3
            end
        end
    end
end
