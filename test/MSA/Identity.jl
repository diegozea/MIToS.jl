@testset "Identity" begin

    @testset "Float64" begin

        @test_throws ErrorException percentidentity(res"AH", res"AGH")

        @test percentidentity(res"AH", res"AH") == 100.
        @test percentidentity(res"AH", res"AG") == 50.
        @test percentidentity(res"AH", res"RG") == 0.
        @test percentidentity(res"AH-", res"AG-") == 50.
        @test percentidentity(res"A--", res"AG-") == 50.
        # Columns with XAA aren't counted
        @test percentidentity(res"AXA-", res"AG--") == 100 .* (1+0+0+0) / (1+0+1+0)
        @test percentidentity(res"AAX-", res"AG--") == 100 .* (1+0+0+0) / (1+1+0+0)
        @test percentidentity(res"AH-", res"AX-") == 100.
        @test percentidentity(res"AH-", res"XG-") == 0.
        @test percentidentity(res"AGG", res"AHX") == 50.

        @test isnan(percentidentity(res"---", res"---"))
        @test isnan(percentidentity(res"XX-", res"-XX"))
        @test isnan(percentidentity(res"XXX", res"XXX"))
    end

    @testset "MSA 2x3" begin
        αβ = (1,20,21,22)
        for i₁ in αβ, j₁ in αβ, k₁ in αβ
            seq₁ = Residue[i₁,j₁,k₁]
            for i₂ in αβ, j₂ in αβ, k₂ in αβ
                seq₂ = Residue[i₂,j₂,k₂]
                val = percentidentity(seq₁, seq₂)
                if  sum((seq₁ .== GAP) .& (seq₂ .== GAP)) +
                    sum((seq₁ .== XAA) .| (seq₂ .== XAA)) == 3

                    @test isnan(val)
                elseif seq₁ == seq₂
                    @test val == 100.0
                else
                    @test 100.0 >= val >= 0.0
                end
            end
        end
    end

    @testset "Bool" begin

        @test percentidentity(res"A--", res"AG-", 40.)
        @test !percentidentity(res"A--", res"AG-",60)
        # Columns with XAA aren't counted
        @test percentidentity(res"AXA-", res"AG--", 40.)
        @test !percentidentity(res"AAX-", res"AG--", 60.)
        @test percentidentity(res"AH-", res"AX-", 100.)
        @test percentidentity(res"AH-", res"XG-", 0.)
        @test !percentidentity(res"AGG", res"AHX", 62.) # 50% < 62%
    end

    @testset "MSA" begin

        fasta = read(joinpath(pwd(), "data", "Gaoetal2011.fasta"), FASTA)
        id = percentidentity(fasta)

        @test id[1,1] == 100.0

        @test isapprox(id[1,2], 83.33, atol=0.01)
        @test isapprox(id[1,3], 83.33, atol=0.01)
        @test isapprox(id[3,1], 83.33, atol=0.01)
        @test isapprox(id[1,6], 33.33, atol=0.01)
        @test isapprox(id[4,5], 83.33, atol=0.01)

        @test id[5,6] == 100.0

        @test maximum(id) == 100.0
        @test isapprox(minimum(id), 33.33, atol=0.01)

        @testset "SequenceIdentityMatrix" begin

            @test isapprox(sum(id[:,3]), 100.0 + 50.0 + 2*(200/6) + 2*(500/6))
            id[4,3] = 80
            @test isapprox(sum(id[:,3]), 100.0 + 80.0 + 2*(200/6) + 2*(500/6))
        end

        @testset "Gaps" begin

            aln = read(joinpath(pwd(), "data", "gaps.txt"), Raw)
            id = percentidentity(aln)

            @test id[1,1] == 100.0
            @test id[1,2] == 90.0
            @test id[1,3] == 80.0
            @test id[2,3] ≈ 800/9
        end

        @testset "Mean percent identity" begin

            aln = Residue[ '-' '-' 'G' 'G' 'G' '-'
                           '-' '-' '-' 'G' 'G' 'G' ]
            # identities    0   0   0   1   1   0 sum 2
            # aligned res   0   0   1   1   1   1 sum 4
            @test percentidentity(aln)[1, 2] == 50.0 # 2 / 4
            @test meanpercentidentity(aln)   == 50.0

            msa    = rand(Residue, 400, 2)
            msa300 = msa[1:300, :]

            @test mean(getlist(percentidentity(msa300))) == meanpercentidentity(msa300)
            @test isapprox( mean(getlist(percentidentity(msa))),
                            meanpercentidentity(msa), atol=0.5 )
            @test mean(getlist(percentidentity(msa))) == meanpercentidentity(msa,exact=true)
        end
    end

    @testset "Percent Similarity" begin

        fasta = read(joinpath(pwd(), "data", "Gaoetal2011.fasta"), FASTA)

        @testset "Gaps" begin

            @test percentsimilarity(res"AH", res"IM") == 50.0

            @test percentsimilarity(res"-AH", res"--H") == 50.0
            @test percentsimilarity(res"-AH", res"-AH") == 100.0
            @test percentsimilarity(res"-AH",res"-AH") ==
                percentsimilarity(res"AH", res"AH")
            @test percentsimilarity(res"-AH",res"-AH") ==
                percentsimilarity(res"--AH", res"--AH")
            # Residues outside the alphabet aren't counted (i.e.: XAA)
            @test percentsimilarity(res"-AH",res"-AH") ==
                percentsimilarity(res"XXAH", res"AAAH")
                @test percentsimilarity(res"AH",res"AH") ==
                    percentsimilarity(res"GGAH", res"AAAH",
                        ReducedAlphabet("AH"))
        end

        @testset "Using SMS's Ident and Sim residue groups" begin

            sim = percentsimilarity(fasta,
                  ReducedAlphabet("(GAVLI)(FYW)(ST)(KRH)(DENQ)P(CM)"))

            @test eltype(sim) == Float64
            @test sim[1,1] == 100.0

            @test isapprox(sim[1,2],  83.33, atol=0.01)
            @test isapprox(sim[1,3], 100.00, atol=0.01)
            @test isapprox(sim[1,4],  66.67, atol=0.01)
            @test isapprox(sim[1,5],  50.00, atol=0.01)
            @test isapprox(sim[1,6],  50.00, atol=0.01)
            @test isapprox(sim[2,3],  83.33, atol=0.01)
            @test isapprox(sim[2,4],  50.00, atol=0.01)
            @test isapprox(sim[2,5],  50.00, atol=0.01)
            @test isapprox(sim[2,6],  50.00, atol=0.01)
            @test isapprox(sim[3,4],  66.67, atol=0.01)
            @test isapprox(sim[3,5],  50.00, atol=0.01)
            @test isapprox(sim[3,6],  50.00, atol=0.01)
            @test isapprox(sim[4,5],  83.33, atol=0.01)
            @test isapprox(sim[4,6],  83.33, atol=0.01)
            @test isapprox(sim[5,6], 100.00, atol=0.01)
        end

        @testset "Using Bio3D's (2.2) seqidentity residue groups" begin

            sim = percentsimilarity(fasta,
                  ReducedAlphabet("(GA)(MVLI)(FYW)(ST)(KRH)(DE)(NQ)PC"), out=Float16)
            bio3d = [   1.000 0.833 1.000 0.667 0.500 0.500
                        0.833 1.000 0.833 0.500 0.500 0.500
                        1.000 0.833 1.000 0.667 0.500 0.500
                        0.667 0.500 0.667 1.000 0.833 0.833
                        0.500 0.500 0.500 0.833 1.000 1.000
                        0.500 0.500 0.500 0.833 1.000 1.000 ] .* 100.00

            @test eltype(sim) == Float16

            for i in 1:6
                for j in 1:6
                    @test isapprox(sim[i,j], bio3d[i,j], atol=0.1)
                end
            end
        end

        @testset "Percent identity" begin

            for αβ in (UngappedAlphabet(), GappedAlphabet())

                @test_throws ErrorException percentsimilarity(res"AH", res"AGH", αβ)

                @test percentsimilarity(res"AH", res"AH", αβ) == 100.
                @test percentsimilarity(res"AH", res"AG", αβ) == 50.
                @test percentsimilarity(res"AH", res"RG", αβ) == 0.
                @test percentsimilarity(res"AH-",res"AG-", αβ) == 50.
                @test percentsimilarity(res"A--",res"AG-", αβ) == 50.
                # Columns with XAA aren't counted
                @test percentsimilarity(res"AXA-",res"AG--", αβ) == 100 .* (1+0+0+0) / (1+0+1+0)
                @test percentsimilarity(res"AAX-",res"AG--", αβ) == 100 .* (1+0+0+0) / (1+1+0+0)
                @test percentsimilarity(res"AH-", res"AX-", αβ) == 100.
                @test percentsimilarity(res"AH-", res"XG-", αβ) == 0.
                @test percentsimilarity(res"AGG", res"AHX", αβ) == 50.

                @test isnan(percentsimilarity(res"---", res"---", αβ))
                @test isnan(percentsimilarity(res"XX-", res"-XX", αβ))
                @test isnan(percentsimilarity(res"XXX", res"XXX", αβ))
            end
        end
    end
end
