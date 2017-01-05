@testset "Pairwise Gap Percentage" begin

    @testset "Simple" begin
        file = joinpath(pwd(), "data", "simple.fasta")
        mat = [ 0. 0.
                0. 0. ]

        (gu, gi) = pairwisegapfraction(file, FASTA)
        @test gu == mat
        @test gi == mat
    end

    @testset "Gaps" begin

        file = joinpath(pwd(), "data", "gaps.txt")
        cl = hobohmI(read(file, Raw), 62)
        gu, gi = pairwisegapfraction(file, Raw)
        ncl = nclusters(cl)

        @test gu[1, 1] ≈ 0.0
        @test gi[1, 1] ≈ 0.0
        @test gu[1, 2] ≈ 100.0 * getweight(cl, 10)/ncl
        @test gi[1, 2] ≈ 0.0
        @test gu[10, 9] ≈ 100.0 * (ncl - getweight(cl, 1))/ncl
        @test gi[10, 9] ≈ 100.0 * (ncl - getweight(cl, 1) - getweight(cl, 2))/ncl
        @test gu[10, 10] ≈ 100.0 * (ncl - getweight(cl, 1))/ncl
        @test gu[10, 10] ≈ 100.0 * (ncl - getweight(cl, 1))/ncl
    end
end
