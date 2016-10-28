@testset "Identity" begin

    @testset "Float64" begin

        @test_throws ErrorException percentidentity(res"AH", res"AGH")

        @test percentidentity(res"AH", res"AH") == 100.
        @test percentidentity(res"AH", res"AG") == 50.
        @test percentidentity(res"AH", res"RG") == 0.
        @test percentidentity(res"AH-", res"AG-") == 50.
        @test percentidentity(res"A--", res"AG-") == 50.
        # Columns with XAA aren't counted
        @test percentidentity(res"AXA-", res"AG--") == 100.*(1+0+0+0)/(1+0+1+0)
        @test percentidentity(res"AAX-", res"AG--") == 100.*(1+0+0+0)/(1+1+0+0)
        @test percentidentity(res"AH-", res"AX-") == 100.
        @test percentidentity(res"AH-", res"XG-") == 0.
        @test percentidentity(res"AGG", res"AHX") == 50.
    end

    @testset "Bool" begin

        @test percentidentity(res"A--", res"AG-", 40.)
        @test !percentidentity(res"A--", res"AG-",60)
        # Columns with XAA aren't counted
        @test percentidentity(res"AXA-", res"AG--", 40.)
        @test !percentidentity(res"AAX-", res"AG--", 60.)
        @test percentidentity(res"AH-", res"AX-", 100.)
        @test percentidentity(res"AH-", res"XG-", 0.)
        @test percentidentity(res"AGG", res"AHX", 62.)
    end

    @testset "MSA" begin

        fasta = read(joinpath(pwd(), "data", "Gaoetal2011.fasta"), FASTA)
        id = percentidentity(fasta)

        @test id[1,1] == 100.0

        @test id[1,2] ≈ 83.33
        @test id[1,3] ≈ 83.33
        @test id[3,1] ≈ 83.33
        @test id[1,6] ≈ 33.33
        @test id[4,5] ≈ 83.33

        @test id[5,6] == 100.0

        @test maximum(id) == 100.0
        @test minimum(id) ≈ 33.33

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

    end
end
