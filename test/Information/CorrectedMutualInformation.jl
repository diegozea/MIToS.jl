@testset "CorrectedMutualInformation" begin

    Gaoetal2011 = joinpath(pwd(), "data", "Gaoetal2011.fasta")

    function gao11_buslje09(measure)
        filename = string("data_Gaoetal2011_soft_Busljeetal2009_measure_", measure, ".txt")
        joinpath(pwd(), "data", filename)
    end

    ## Column numbers for the output of Buslje et. al. 2009
    SCORE = 9
    ZSCORE = 12

    MIToS_SCORE = 2
    MIToS_ZSCORE = 1

    @testset "APC!" begin

        for MI in ([ 0. 2. 4.
                     2. 0. 6.
                     4. 6. 0. ],
                    PairwiseListMatrix([2., 4., 6.]),
                    PairwiseListMatrix([0., 2., 4., 0., 6., 0.], true),
                    NamedArray(PairwiseListMatrix([2., 4., 6.])),
                    NamedArray(PairwiseListMatrix([0., 2., 4., 0., 6., 0.],true)))

            if isa(MI, Matrix{Float64})
                @test Information._mean_column(MI) == [3., 4., 5.]
                @test Information._mean_total(MI) == 4.
            end

            MIp = APC!(MI)

            @test filter(x -> !isnan(x), vec(Matrix(MIp))) ≈ filter(x -> !isnan(x), vec(
                [   NaN -1.0  0.25
                    -1.0  NaN  1.00
                    0.25 1.00   NaN ]))
        end
    end

    @testset "Simple example I" begin

        aln = Residue[ 'A' 'A'
                       'A' 'R' ]

        @testset "MI" begin

            # Fill the Pij matrix of the example
            Pij = zeros(Float64, 20, 20);
            N = 2 # There are 2 sequences
            Pij[1,1] = (1/N) # A A
            Pij[1,2] = (1/N) # A R
            @test sum(Pij) ≈ 1.0 # Is the matrix correct?
            # Fill the marginals
            Pi = dropdims(Base.sum(Pij, dims=2), dims=2)
            @test Pi[1] == 1.0 # Is always A
            Pj = dropdims(Base.sum(Pij, dims=1), dims=1)
            @test Pj[1] == 0.5 # A
            @test Pj[2] == 0.5 # R
            # Start to sum with 0.0
            total = 0.0
            for i in 1:20, j in 1:20
                if Pij[i,j] != 0.0 && Pi[i] != 0.0 && Pj[j] != 0.0
                    total += (Pij[i,j] * log(Pij[i,j]/(Pi[i]*Pj[j]))) # 0.5 * log(0.5/(0.5 * 1.0)) == 0.0
                end
            end
            # Compare with MIToS result
            mi = mapcolpairfreq!(mutual_information, aln,
                    Probabilities(ContingencyTable(Float64,Val{2},UngappedAlphabet())))
            @test mi[1,2] == total
        end

        @testset "MI: Using pseudocount (0.05)" begin

            Pij = zeros(Float64, 20, 20);
            N = (400 * 0.05) + 2 # 0.05 is the pseudocount and there are 2 sequences
            fill!(Pij, 0.05/N);
            Pij[1,1] =(1.05/N) # A A
            Pij[1,2] =(1.05/N) # A R
            @test sum(Pij) ≈ 1.0
            Pi = dropdims(Base.sum(Pij, dims=2), dims=2)
            Pj = dropdims(Base.sum(Pij, dims=1), dims=1)
            total = 0.0
            for i in 1:20, j in 1:20
                total += (Pij[i,j] * log(Pij[i,j]/(Pi[i]*Pj[j])))
            end
            # Compare with MIToS result
            mi = mapcolpairfreq!(mutual_information, aln,
                    Probabilities(ContingencyTable(Float64,Val{2},UngappedAlphabet())),
                    pseudocounts=AdditiveSmoothing(0.05))
            @test mi[1,2] ≈ total

            @testset "APC!" begin

                APC!(mi)
                @test isnan(mi[1,1])
                @test isapprox(mi[1,2], 0.0, rtol=1e-16)

                zerodiagonal = Float64[ 0.0 total
                                        total 0.0 ]

                # MI.. == MI.j == MIi. == total
                # APC = (total * total) / total == total
                # MI - APC == total - total == 0.0
                APC!(zerodiagonal)
                @test zerodiagonal[1,2] ≈ 0.0
            end
        end

        @testset "Z-score" begin

            # 2 Possibilities:  A A   A A
            #                   A R   R A
            # Is almost the same MI, Z score should be 0.0
            results = buslje09(aln, lambda=0.05, clustering=false, apc=false)

            @test results[MIToS_ZSCORE][1,2] ≈ 0.0

            @test aln == Residue[ 'A' 'A'
                                  'A' 'R' ]
        end
    end

    @testset "Simple example II" begin

        aln = Residue[ 'R' 'A'
                       'A' 'R' ]

        @testset "MI" begin

            Pij = zeros(Float64, 20, 20);
            N = 2 # There are 2 sequences
            Pij[2,1] = (1/N) # R A
            Pij[1,2] = (1/N) # A R
            @test sum(Pij) ≈ 1.0
            Pi = dropdims(Base.sum(Pij, dims=2), dims=2)
            @test Pi[2] == 0.5 # R
            @test Pi[1] == 0.5 # A
            Pj = dropdims(Base.sum(Pij, dims=1), dims=1)
            @test Pj[1] == 0.5 # A
            @test Pj[2] == 0.5 # R
            total = 0.0
            for i in 1:20, j in 1:20
                if Pij[i,j] != 0.0 && Pi[i] != 0.0 && Pj[j] != 0.0
                    total += (Pij[i,j] * log(Pij[i,j]/(Pi[i]*Pj[j])))
                end
            end
            @test total == 0.5 * log(2) + 0.5 * log(2) # 0.5/(0.5*0.5) == 2
            # Compare with MIToS result
            mi = mapcolpairfreq!(mutual_information, aln,
                    Probabilities(ContingencyTable(Float64,Val{2},UngappedAlphabet())))
            @test mi[1,2] == total
        end

        @testset "MI: Using pseudocount (0.05)" begin

            Pij = zeros(Float64, 20, 20);
            N = (400 * 0.05) + 2 # 0.05 is the pseudocount and there are 2 sequences
            fill!(Pij, 0.05/N);
            Pij[2,1] =(1.05/N) # R A
            Pij[1,2] =(1.05/N) # A R
            @test sum(Pij) ≈ 1.0
            Pi = dropdims(Base.sum(Pij, dims=2), dims=2)
            Pj = dropdims(Base.sum(Pij, dims=1), dims=1)
            total = 0.0
            for i in 1:20, j in 1:20
                total += (Pij[i,j] * log(Pij[i,j]/(Pi[i]*Pj[j])))
            end
            # Compare with MIToS result
            mi = mapcolpairfreq!(mutual_information, aln,
                    Probabilities(ContingencyTable(Float64,Val{2},UngappedAlphabet())),
                    pseudocounts=AdditiveSmoothing(0.05))
            @test mi[1,2] ≈ total

            @testset "APC!" begin

                APC!(mi)
                @test isnan(mi[1,1])
                @test mi[1,2] ≈ 0.0

                zerodiagonal = Float64[ 0.0 total
                                        total 0.0 ]

                # MI.. == MI.j == MIi. == total
                # APC = (total * total) / total == total
                # MI - APC == total - total == 0.0
                APC!(zerodiagonal)
                @test zerodiagonal[1,2] ≈ 0.0
            end
        end

        @testset "Z-score" begin

            # 4 Possibilities:  R A   R A   A R   A R
            #                   A R   R A   A R   R A
            mi = mapcolpairfreq!(mutual_information, aln,
                    Probabilities(ContingencyTable(Float64,Val{2},UngappedAlphabet())),
                    pseudocounts=AdditiveSmoothing(0.05))
            other_posib = Residue[ 'A' 'R'
                                   'A' 'R' ]
            other_mi = mapcolpairfreq!(mutual_information, other_posib,
                        Probabilities(ContingencyTable(Float64,Val{2},UngappedAlphabet())),
                        pseudocounts=AdditiveSmoothing(0.05))
            # The mean should be:
            r_mean = 0.5 * ( mi[1,2] + other_mi[1,2] )
            # And the std should be similar to:
            r_std = 0.5 * ( sqrt((mi[1,2]-r_mean)^2) + sqrt((other_mi[1,2]-r_mean)^2) )

            results = buslje09(aln, lambda=0.05, clustering=false, apc=false, samples=100)
            @test isapprox(results[MIToS_ZSCORE][1,2], (mi[1,2]-r_mean)/r_std, atol=0.55)

            results = buslje09(aln, lambda=0.05, clustering=false, apc=false, samples=1000)
            @test isapprox(results[MIToS_ZSCORE][1,2], (mi[1,2]-r_mean)/r_std, atol=0.15)

            results = buslje09(aln, lambda=0.05, clustering=false, apc=false, samples=10000)
            @test isapprox(results[MIToS_ZSCORE][1,2], (mi[1,2]-r_mean)/r_std, atol=0.055)

            @test aln == Residue[ 'R' 'A'
                                  'A' 'R' ]
        end
    end

    @testset "Results from Buslje et. al. 2009" begin

        @testset "Simple" begin
            data = readdlm(joinpath(pwd(), "data",
                           "data_simple_soft_Busljeetal2009_measure_MI.txt"), comments=true)
            results = buslje09(joinpath(pwd(), "data", "simple.fasta"), FASTA,
                               lambda=0.0, clustering=false, apc=false)

            @test isapprox(Float64(data[1, SCORE]), results[MIToS_SCORE][1,2], atol=1e-6)
            @test isapprox(Float64(data[1, ZSCORE]), results[MIToS_ZSCORE][1,2], atol=2.)
        end

        @testset "Gaoetal2011" begin
            data = readdlm(gao11_buslje09("MI"), comments=true)
            results = buslje09(Gaoetal2011, FASTA, lambda=0.0, clustering=false, apc=false)

            @test isapprox([ Float64(x) for x in data[:,SCORE] ],
                           matrix2list(results[MIToS_SCORE]), atol=1e-6)
            @test isapprox([ Float64(x) for x in data[:,ZSCORE] ],
                           matrix2list(results[MIToS_ZSCORE]), atol=2.)

            # println(cor(convert(Vector{Float64},data[:,ZSCORE]),matrix2list(results[MIToS_ZSCORE])))

            @testset "cMI" begin
                result = copy(results[1])
                result[diagind(result)] .= 0.0
                @test cumulative(results[1], -Inf) == sum(result, dims=1)
                @test all((cumulative(results[2],0.5) .-
                           [0.0,0.0,1.38629,1.38629,1.38629,0.0]') .< 0.00001)
                           #        1.38629 = 0.693147 + 0.693147
            end

            result_0_05 = buslje09(Gaoetal2011,FASTA,lambda=0.05,clustering=false,apc=false)
            @test isapprox(result_0_05[MIToS_SCORE][1,2], 0.33051006116310444, atol=1e-14)
        end
    end

    @testset "MI + clustering" begin

        data = readdlm(gao11_buslje09("MI_clustering"), comments=true)
        results = buslje09(Gaoetal2011, FASTA, lambda=0.0, clustering=true, apc=false)

        @test isapprox([ Float64(x) for x in data[:,SCORE] ],
                       matrix2list(results[MIToS_SCORE]), atol=1e-6)
        @test isapprox([ Float64(x) for x in data[:,ZSCORE] ],
                       matrix2list(results[MIToS_ZSCORE]), atol=2.)
    end

    @testset "MIp" begin

        data = readdlm(gao11_buslje09("MI_APC"), comments=true)
        results = buslje09(Gaoetal2011, FASTA, lambda=0.0, clustering=false, apc=true)

        @test isapprox([ Float64(x) for x in data[:,SCORE] ],
                       matrix2list(results[MIToS_SCORE]), atol=1e-6)
        @test isapprox([ Float64(x) for x in data[:,ZSCORE] ],
                       matrix2list(results[MIToS_ZSCORE]), atol=2.0)

        @test isapprox(results[MIToS_SCORE][5,6], 0.018484, atol=0.000001)
    end

    @testset "MIp + clustering" begin

        data = readdlm(gao11_buslje09("MI_APC_clustering"), comments=true)
        results = buslje09(Gaoetal2011, FASTA, lambda=0.0, clustering=true, apc=true)

        @test isapprox([ Float64(x) for x in data[:,SCORE] ],
                       matrix2list(results[MIToS_SCORE]), atol=1e-6)
        @test isapprox([ Float64(x) for x in data[:,ZSCORE] ],
                       matrix2list(results[MIToS_ZSCORE]), atol=2.)

        @test isapprox(results[MIToS_SCORE][5,6], 0.018484, atol=0.000001)
    end

    @testset "BLMI" begin

        @testset "Simple" begin
            file = joinpath(pwd(), "data", "simple.fasta")
            busl = buslje09(file, FASTA)
            blmi = BLMI(file, FASTA)

            @test PairwiseListMatrices.getlist(busl[1]) ≈
                  PairwiseListMatrices.getlist(blmi[1])
            @test PairwiseListMatrices.getlist(busl[2]) ≈
                  PairwiseListMatrices.getlist(blmi[2])
        end

        @testset "Gaoetal2011" begin
            msa  = read(Gaoetal2011, FASTA)
            busl = buslje09(Gaoetal2011, FASTA, lambda=0.0, samples=0)
            blmi = BLMI(msa, lambda=0.0, beta=0.0, samples=5)
            # BLMI should be equal to Buslje09 if beta is zero

            @test PairwiseListMatrices.getlist(busl[2]) ≈
                  PairwiseListMatrices.getlist(blmi[2]) # MIapc
            @test msa == read(Gaoetal2011, FASTA)
        end

        @testset "Gaoetal2011, lambda 0.05" begin

            busl = buslje09(Gaoetal2011, FASTA, lambda=0.5, samples=0)
            blmi = BLMI(Gaoetal2011, FASTA, lambda=0.5, beta=0.0, samples=5)
            # BLMI should be equal to Buslje09 if beta is zero

            @test PairwiseListMatrices.getlist(busl[2]) ≈
                  PairwiseListMatrices.getlist(blmi[2]) # MIapc
        end
    end
end
