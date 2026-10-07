@testset "Iterations" begin

    function _test_same_named_plm(a, b)
        @test typeof(a) === typeof(b)
        @test size(a) == size(b)
        @test isequal(getlist(getarray(a)), getlist(getarray(b)))
        @test isequal(getdiag(getarray(a)), getdiag(getarray(b)))
        @test dimnames(a) == dimnames(b)
        @test names(a, 1) == names(b, 1)
        @test names(a, 2) == names(b, 2)
    end

    function _test_same_pair_table(a, b)
        @test isequal(gettable(a), gettable(b))
        @test isequal(getmarginals(a), getmarginals(b))
        @test isequal(gettotal(a), gettotal(b))
        @test isequal(getcontingencytable(a).temporal, getcontingencytable(b).temporal)
    end

    function _test_threaded_pairs(f, msa, table; kwargs...)
        serial_table = deepcopy(table)
        threaded_table = deepcopy(table)
        storage = gettablearray(threaded_table)
        serial = mapcolpairfreq!(f, msa, serial_table; kwargs...)
        threaded = mapcolpairfreq!(f, msa, threaded_table; threads = true, kwargs...)
        _test_same_named_plm(serial, threaded)
        _test_same_pair_table(serial_table, threaded_table)
        @test storage === gettablearray(threaded_table)
        threaded
    end

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
            Frequencies(ContingencyTable(Float64, Val{2}, UngappedAlphabet())),
            usediagonal = false,
        )
        nmi_mat = convert(Matrix{Float64}, getarray(nmi))
        @test isapprox(nmi_mat, result, rtol = 1e-4)

        nmi_t = mapseqpairfreq!(
            normalized_mutual_information,
            permutedims(aln),
            Frequencies(ContingencyTable(Float64, Val{2}, UngappedAlphabet())),
            usediagonal = false,
        )
        @test nmi_mat == convert(Matrix{Float64}, getarray(nmi_t))
    end

    @testset "Threaded mapcolpairfreq!" begin
        synthetic = Residue[
            'A' 'A' 'R' 'N'
            'A' 'R' 'R' 'N'
            'R' 'A' 'A' 'D'
            'R' 'R' 'A' 'D'
            'A' 'A' 'A' 'N'
        ]

        synthetic_serial = mapcolpairfreq!(
            Information._mutual_information,
            synthetic,
            Probabilities(ContingencyTable(Float64, Val{2}, UngappedAlphabet()));
            usediagonal = true,
            pseudocounts = AdditiveSmoothing(0.05),
        )
        synthetic_threaded = mapcolpairfreq!(
            Information._mutual_information,
            synthetic,
            Probabilities(ContingencyTable(Float64, Val{2}, UngappedAlphabet()));
            usediagonal = true,
            pseudocounts = AdditiveSmoothing(0.05),
            threads = true,
        )

        _test_same_named_plm(synthetic_serial, synthetic_threaded)

        aln = read_file(joinpath(DATA, "Gaoetal2011.fasta"), FASTA)
        serial = mapcolpairfreq!(
            normalized_mutual_information,
            aln,
            Frequencies(ContingencyTable(Float64, Val{2}, UngappedAlphabet()));
            usediagonal = false,
        )
        threaded = mapcolpairfreq!(
            normalized_mutual_information,
            aln,
            Frequencies(ContingencyTable(Float64, Val{2}, UngappedAlphabet()));
            usediagonal = false,
            threads = true,
        )

        _test_same_named_plm(serial, threaded)
    end

    @testset "Threaded pair options and table state" begin
        rng = Random.MersenneTwister(197)
        residues = rand(rng, res"ARNDCQEGHILKMFPSTWYV-X", 8, 9)
        residues[:, 1] .= GAP
        msa = NamedArray(
            residues,
            (["seq_$i" for i = 1:8], ["site_$(2i)" for i = 1:9]),
            ("Seq", "Col"),
        )
        weighting = (
            NoClustering(),
            Weights([0.1, 0.3, 0.7, 1.0, 0.2, 0.6, 0.4, 0.9]),
            hobohmI(residues, 62),
        )
        for wrapper in (Frequencies, Probabilities),
            alphabet in (
                UngappedAlphabet(),
                GappedAlphabet(),
                ReducedAlphabet("(AILMV)(RHK)(NQST)(DE)(FWY)CGP"),
            ),
            weights in weighting,
            usediagonal in (false, true)

            table = wrapper(ContingencyTable(Float64, Val{2}, alphabet))
            fill!(getcontingencytable(table), 2.0)
            _test_threaded_pairs(
                mutual_information,
                msa,
                table;
                weights = weights,
                pseudocounts = AdditiveSmoothing(0.05),
                usediagonal = usediagonal,
                diagonalvalue = NaN,
                base = 2,
            )
        end

        for usediagonal in (false, true)
            _test_threaded_pairs(
                mutual_information,
                residues,
                Probabilities(ContingencyTable(Float64, Val{2}, UngappedAlphabet()));
                pseudocounts = AdditiveSmoothing(0.05),
                pseudofrequencies = BLOSUM_Pseudofrequencies(8.0, 8.512),
                usediagonal = usediagonal,
                base = 2,
            )
            _test_threaded_pairs(
                normalized_mutual_information,
                residues,
                Frequencies(ContingencyTable(Float32, Val{2}, GappedAlphabet()));
                usediagonal = usediagonal,
                diagonalvalue = -1.0f0,
            )
        end
        _test_same_named_plm(
            normalized_mutual_information(msa),
            normalized_mutual_information(msa; threads = true),
        )
    end

    @testset "Threaded pair boundaries and nonfinite scores" begin
        for nseq in (0, 5), ncol in (0, 1, 2, 3, 4, 7, 11), usediagonal in (false, true)
            msa = rand(Random.MersenneTwister(ncol), res"ARND-", nseq, ncol)
            table = Frequencies(ContingencyTable(Float64, Val{2}, UngappedAlphabet()))
            fill!(getcontingencytable(table), 2.0)
            _test_threaded_pairs(
                mutual_information,
                msa,
                table;
                usediagonal = usediagonal,
                diagonalvalue = -0.0,
            )
        end
        for residue in (GAP, Residue('A'))
            _test_threaded_pairs(
                normalized_mutual_information,
                fill(residue, 5, 7),
                Probabilities(ContingencyTable(Float64, Val{2}, UngappedAlphabet())),
            )
        end
    end

    @testset "Yielding and nested threaded callbacks" begin
        msa = rand(Random.MersenneTwister(197), res"ARNDCQEGHILKMFPSTWYV-", 40, 20)
        table = Frequencies(ContingencyTable(Float64, Val{2}, GappedAlphabet()))
        function yielding_score(t; base = 2)
            before = copy(gettablearray(t))
            for _ = 1:8
                yield()
            end
            isequal(before, gettablearray(t)) || error("Scratch table shared across tasks")
            mutual_information(t; base = base)
        end
        for usediagonal in (false, true)
            serial = mapcolpairfreq!(
                mutual_information,
                msa,
                deepcopy(table);
                usediagonal = usediagonal,
                base = 2,
            )
            # Concurrent calls from worker tasks also exercise nested parallelism.
            tasks = [
                Threads.@spawn(
                    mapcolpairfreq!(
                        yielding_score,
                        msa,
                        deepcopy(table);
                        usediagonal = usediagonal,
                        threads = true,
                        base = 2,
                    )
                ) for _ = 1:8
            ]
            for task in tasks
                _test_same_named_plm(serial, fetch(task))
            end
        end

        # Ensure a multi-worker test run actually exercises parallel callbacks.
        owners = Set{Task}()
        owner_lock = ReentrantLock()
        function record_task(t)
            lock(owner_lock) do
                push!(owners, current_task())
            end
            gettotal(t)
        end
        mapcolpairfreq!(record_task, msa, table; threads = true)
        @test length(owners) == min(Threads.nthreads(), 210)
    end

    @testset "Gaps" begin

        function _gaps(
            table::Union{
                Probabilities{Float64,1,GappedAlphabet},
                Frequencies{Float64,1,GappedAlphabet},
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

        colcount = mapcolfreq!(_gaps, gaps, Frequencies(table))
        @test all((vec(getarray(colcount)) .- ngaps) .== 0.0)

        colfract = mapcolfreq!(_gaps, gaps, Probabilities(table))
        @test all((vec(getarray(colfract)) .- ngaps ./ 10.0) .== 0.0)

        seqcount = mapseqfreq!(_gaps, gaps, Frequencies(table))
        @test all((vec(getarray(seqcount)) .- ngaps) .== 0.0)

        seqfract = mapseqfreq!(_gaps, gaps, Probabilities(table))
        @test all((vec(getarray(seqfract)) .- ngaps ./ 10.0) .== 0.0)
    end

    @testset "Passing keyword arguments" begin

        # Dummy function that returns the value of the keyword argument `karg` for testing
        f(table; karg::Float64 = 0.0) = karg

        msa = rand(Random.MersenneTwister(123), Residue, 4, 10)
        table_1d = Frequencies(ContingencyTable(Float64, Val{1}, UngappedAlphabet()))
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
