@testset "Iterations" begin

    @testset "NMI" begin
        # This is the example of MI(X, Y)/H(X, Y) from:
        #
        # Gao, H., Dou, Y., Yang, J., & Wang, J. (2011).
        # New methods to measure residues coevolution in proteins.
        # BMC bioinformatics, 12(1), 206.

        aln = read_file(joinpath(DATA, "Gaoetal2011.fasta"), FASTA)
        result = Float64[
            0 0 0 0 0 0
            0 0 0 0 0 0
            0 0 0 1 1 0.296
            0 0 1 0 1 0.296
            0 0 1 1 0 0.296
            0 0 0.296 0.296 0.296 0
        ]

        nmi = mapcolpairfreq!(
            normalized_mutual_information,
            aln,
            Counts(ContingencyTable(Float64, Val{2}, UngappedAlphabet())),
            usediagonal = false,
        )
        nmi_mat = convert(Matrix{Float64}, getarray(nmi))
        @test isapprox(nmi_mat, result, rtol = 1e-4)

        nmi_t = mapseqpairfreq!(
            normalized_mutual_information,
            permutedims(aln),
            Counts(ContingencyTable(Float64, Val{2}, UngappedAlphabet())),
            usediagonal = false,
        )
        @test nmi_mat == convert(Matrix{Float64}, getarray(nmi_t))
    end

    @testset "Gaps" begin

        function _gaps(
            table::Union{
                Probabilities{Float64,1,GappedAlphabet},
                Counts{Float64,1,GappedAlphabet},
            },
        )
            table[GAP]
        end

        table = ContingencyTable(Float64, Val{1}, GappedAlphabet())

        gaps = read_file(joinpath(DATA, "gaps.txt"), Raw)

        # THAYQAIHQV 0
        # THAYQAIHQ- 0.1
        # THAYQAIH-- 0.2
        # THAYQAI--- 0.3
        # THAYQA---- 0.4
        # THAYQ----- 0.5
        # THAY------ 0.6
        # THA------- 0.7
        # TH-------- 0.8
        # T--------- 0.9

        ngaps = Float64[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

        colcount = mapcolfreq!(_gaps, gaps, Counts(table))
        @test all((vec(getarray(colcount)) .- ngaps) .== 0.0)

        colfract = mapcolfreq!(_gaps, gaps, Probabilities(table))
        @test all((vec(getarray(colfract)) .- ngaps ./ 10.0) .== 0.0)

        seqcount = mapseqfreq!(_gaps, gaps, Counts(table))
        @test all((vec(getarray(seqcount)) .- ngaps) .== 0.0)

        seqfract = mapseqfreq!(_gaps, gaps, Probabilities(table))
        @test all((vec(getarray(seqfract)) .- ngaps ./ 10.0) .== 0.0)
    end

    @testset "Passing keyword arguments" begin

        # Dummy function that returns the value of the keyword argument `karg` for testing
        f(table; karg::Float64 = 0.0) = karg

        msa = rand(Random.MersenneTwister(123), Residue, 4, 10)
        table_1d = Counts(ContingencyTable(Float64, Val{1}, UngappedAlphabet()))
        table_2d = Probabilities(ContingencyTable(Float64, Val{2}, UngappedAlphabet()))

        @test sum(mapseqfreq!(f, msa, deepcopy(table_1d))) == 0.0
        @test sum(mapseqfreq!(f, msa, deepcopy(table_1d), karg = 1.0)) == 4.0

        @test sum(mapcolfreq!(f, msa, deepcopy(table_1d))) == 0.0
        @test sum(mapcolfreq!(f, msa, deepcopy(table_1d), karg = 1.0)) == 10.0

        @test sum(mapseqpairfreq!(f, msa, deepcopy(table_2d))) == 0
        @test sum(mapseqpairfreq!(f, msa, deepcopy(table_2d), karg = 1.0)) == 16.0

        @test sum(mapcolpairfreq!(f, msa, deepcopy(table_2d))) == 0
        @test sum(mapcolpairfreq!(f, msa, deepcopy(table_2d), karg = 1.0)) == 100.0
    end

    @testset "mapfreq" begin
        # test mapfreq using the sum function
        msa = rand(Random.MersenneTwister(123), Residue, 4, 10)

        sum_11 = mapfreq(sum, msa, rank = 1, dims = 1)
        @test sum(sum_11) ≈ 4.0
        @test size(sum_11) == (4, 1)

        sum_12 = mapfreq(sum, msa, rank = 1, dims = 2)
        @test sum(sum_12) ≈ 10.0
        @test size(sum_12) == (1, 10)

        sum_21 = mapfreq(sum, msa, rank = 2, dims = 1)
        @test isnan(sum(sum_21))
        @test sum(i for i in sum_21 if !isnan(i)) ≈ 12 # 4 * (4-1)
        @test size(sum_21) == (4, 4)

        sum_22 = mapfreq(sum, msa, rank = 2, dims = 2)
        @test isnan(sum(sum_22))
        @test sum(i for i in sum_22 if !isnan(i)) ≈ 90 # 10 * (10-1)
        @test size(sum_22) == (10, 10)


        # the default is rank = 1 and dims = 2
        @test mapfreq(sum, msa, rank = 1, dims = 2) == mapfreq(sum, msa)

        # probabilities=false
        @test sum(mapfreq(sum, msa, dims = 1, probabilities = false)) ≈ 4 * 10.0
    end
end
