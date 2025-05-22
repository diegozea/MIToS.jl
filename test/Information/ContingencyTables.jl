@testset "ContingencyTables" begin

    @testset "Creation and Getters" begin

        for alphabet in (
            UngappedAlphabet(),
            GappedAlphabet(),
            ReducedAlphabet("(AILMV)(NQST)(RHK)(DE)(FWY)CGP"),
        )
            for N = 1:3
                table = ContingencyTable(Float64, Val{N}, alphabet) # zeros in MIToS 1.0

                @test size(table) == (Int[length(alphabet) for i = 1:N]...,)
                @test length(table) == length(alphabet)^N
                @test length(getmarginals(table)) == length(alphabet) * N
                @test size(getmarginals(table)) == (length(alphabet), N)
                @test sum(gettable(table)) == 0.0
                @test sum(table) == 0.0 # == gettotal(table)
                @test sum(table.temporal) == 0.0
                @test sum(getmarginals(table)) == 0.0
                @test gettotal(table) == 0.0
                @test collect(table) == gettablearray(table) # Iteration interface in MIToS 1.0
                @test isa(gettablearray(table), Array{Float64,N})
                @test isa(getmarginalsarray(table), Array{Float64,2})
            end
        end

        @testset "Similar" begin
            for alphabet in (
                UngappedAlphabet(),
                GappedAlphabet(),
                ReducedAlphabet("(AILMV)(NQST)(RHK)(DE)(FWY)CGP"),
            )
                for N = 1:3
                    table = ContingencyTable(Float64, Val{N}, alphabet) # zeros in MIToS 1.0

                    @test table == similar(table)
                    @test typeof(table) == typeof(similar(table))
                    @test table == similar(table, BigFloat)
                    @test typeof(table) != typeof(similar(table, BigFloat))
                    @test BigFloat == eltype(similar(table, BigFloat))
                end
            end
        end
    end

    @testset "ContingencyTable Show Methods" begin
        # This testset covers show(io, ::MIME"text/plain", ...) for ContingencyTable,
        # Probabilities, and Frequencies.

        alpha = UngappedAlphabet()

        @testset "ContingencyTable show" begin
            # Test N=1
            ct1 = ContingencyTable(Float64, Val{1}, alpha)
            ct1[Residue('A')] = 5.0
            ct1[Residue('R')] = 3.0
            # update_marginals! for N=1 just updates the total based on the table.
            # The `marginals` field itself is not used in the same way as N>1.
            # The `show` method for N=1 does not print a "marginals :" section.
            # We need to ensure total is calculated for the show method.
            # The show method calls `gettotal(table)`, which should be up-to-date if table is modified.
            # Let's assume direct modification of `table.table` and `table.total` or use of internal updates if available.
            # For testing `show`, it's best to have a fully consistent table.
            # `update_marginals!(ct)` is the public API to ensure this.
            update_marginals!(ct1) # Ensures total is correct.

            str_out_ct1 = sprint((io, x) -> show(io, MIME"text/plain"(), x), ct1)

            @test occursin("ContingencyTable{Float64, 1, UngappedAlphabet} :", str_out_ct1)
            @test occursin("\ntable :", str_out_ct1)
            @test occursin("Named Array", str_out_ct1) # Check if NamedArray show is part of it
            @test occursin("Dim_1", str_out_ct1) # Default NamedArray dimension name
            @test occursin(string(Residue('A')), str_out_ct1) # Check for a residue name
            @test !occursin("\n\nmarginals :", str_out_ct1) # Key: No marginals section for N=1
            @test occursin("\n\ntotal : $(gettotal(ct1))", str_out_ct1)
            @test occursin("\n\ntotal : 8.0", str_out_ct1)


            # Test N=2
            ct2 = ContingencyTable(Float64, Val{2}, alpha)
            ct2[Residue('A'), Residue('C')] = 2.0
            ct2[Residue('G'), Residue('T')] = 4.0
            update_marginals!(ct2)

            str_out_ct2 = sprint((io, x) -> show(io, MIME"text/plain"(), x), ct2)

            @test occursin("ContingencyTable{Float64, 2, UngappedAlphabet} :", str_out_ct2)
            @test occursin("\ntable :", str_out_ct2)
            @test occursin("Named Array", str_out_ct2)
            @test occursin("Dim_1", str_out_ct2)
            @test occursin(string(Residue('A')), str_out_ct2)
            @test occursin("\n\nmarginals :", str_out_ct2) # Marginals section for N=2
            # Check for part of marginals NamedArray output
            @test countlines(IOBuffer(str_out_ct2)) > countlines(IOBuffer(str_out_ct1)) # More content due to marginals
            @test occursin("\n\ntotal : $(gettotal(ct2))", str_out_ct2)
            @test occursin("\n\ntotal : 6.0", str_out_ct2)
        end

        @testset "Probabilities show" begin
            ct_for_prob = ContingencyTable(Float64, Val{2}, alpha)
            ct_for_prob[Residue('A'), Residue('A')] = 10.0
            update_marginals!(ct_for_prob)
            
            # Probabilities constructor normalizes the table.
            prob = Probabilities(ct_for_prob, false) # false to not use pseudocounts for simplicity

            str_out_prob = sprint((io, x) -> show(io, MIME"text/plain"(), x), prob)
            str_out_ct_internal = sprint((io, x) -> show(io, MIME"text/plain"(), x), prob.table) # Get show string of internal CT

            expected_prefix = "Probabilities{Float64, 2, UngappedAlphabet} wrapping a "
            @test startswith(str_out_prob, expected_prefix)
            # The rest of the string should be the show output of the wrapped ContingencyTable
            @test strip(str_out_prob[length(expected_prefix)+1:end]) == strip(str_out_ct_internal)
            
            # Verify that the total of the internal CT (which is now probabilities) is 1.0
            @test occursin("\n\ntotal : 1.0", str_out_prob)
        end

        @testset "Frequencies show" begin
            ct_for_freq = ContingencyTable(Int, Val{1}, alpha) # Frequencies often use Int
            ct_for_freq[Residue('T')] = 5
            ct_for_freq[Residue('G')] = 15
            update_marginals!(ct_for_freq)
            
            # Frequencies constructor just wraps the ContingencyTable
            freq = Frequencies(ct_for_freq)

            str_out_freq = sprint((io, x) -> show(io, MIME"text/plain"(), x), freq)
            str_out_ct_internal_freq = sprint((io, x) -> show(io, MIME"text/plain"(), x), freq.table)

            expected_prefix_freq = "Frequencies{Int64, 1, UngappedAlphabet} wrapping a " # Note Int64
            @test startswith(str_out_freq, expected_prefix_freq)
            @test strip(str_out_freq[length(expected_prefix_freq)+1:end]) == strip(str_out_ct_internal_freq)

            # Verify total from the original counts
            @test occursin("\n\ntotal : $(gettotal(ct_for_freq))", str_out_freq)
            @test occursin("\n\ntotal : 20", str_out_freq)
        end
    end

    @testset "Indexing" begin

        for alphabet in (
            UngappedAlphabet(),
            GappedAlphabet(),
            ReducedAlphabet("(AILMV)(NQST)(RHK)(DE)(FWY)CGP"),
        )

            table = ContingencyTable(Float64, Val{2}, alphabet)

            if isa(getalphabet(table), ReducedAlphabet)
                @test table[1, 3] == 0.0
                table[1, 3] = 100.0
                i = 2 * length(getalphabet(table)) + 1 # using one index
                @test table[i] == 100.0
                table[i] = 10.0
                @test table[Residue('A'), Residue('R')] == 10.0 # using two indices
                table[Residue('A'), Residue('R')] = 20.0
                @test table["AILMV", "RHK"] == 20.0
                table["AILMV", "RHK"] = 30.0
                @test table[1, 3] == 30.0
            else
                @test table[1, 2] == 0.0
                table[1, 2] = 100.0
                i = length(getalphabet(table)) + 1 # using one index
                @test table[i] == 100.0
                table[i] = 10.0
                @test table[Residue('A'), Residue('R')] == 10.0 # using two indices
                table[Residue('A'), Residue('R')] = 20.0
                @test table["A", "R"] == 20.0
                table["A", "R"] = 30.0
                @test table[1, 2] == 30.0
            end
        end
    end

    @testset "Update" begin

        for alphabet in (
            UngappedAlphabet(),
            GappedAlphabet(),
            ReducedAlphabet("(AILMV)(NQST)(RHK)(DE)(FWY)CGP"),
        )
            for N = 1:3

                table = ContingencyTable(Float64, Val{N}, alphabet)

                fill!(table.temporal, 1.0)
                @test sum(table.temporal) == 22.0^N
                @test sum(table) == 0.0
                @test sum(getmarginals(table)) == 0.0
                @test gettotal(table) == 0.0

                Information._update!(table)
                @test sum(table.temporal) == 22.0^N
                if isa(getalphabet(table), ReducedAlphabet)
                    @test table[1] == 5.0^N
                    @test sum(table) == 20.0^N
                    @test sum(getmarginals(table)) == N * (20.0^N)
                    @test gettotal(table) == 20.0^N
                    if N == 1
                        @test vec(getarray(getmarginals(table))) ==
                              vec(getarray(gettable(table)))
                    elseif N == 2
                        @test getmarginals(table)[1] ==
                              sum(5.0 * n for n in [5, 4, 3, 2, 3, 1, 1, 1])
                    elseif N == 2
                        @test getmarginals(table)[1] ==
                              sum(5.0 * n * 20.0 for n in [5, 4, 3, 2, 3, 1, 1, 1])
                    end
                else
                    @test table[1] == 1.0
                    @test getmarginals(table)[1] == length(alphabet)^(N - 1)
                    @test sum(table) == length(alphabet)^N
                    @test sum(getmarginals(table)) == float(N) * (length(alphabet)^N)
                    @test gettotal(table) == float(length(alphabet))^N
                end
            end
        end
    end

    @testset "Fill" begin

        for alphabet in (
            UngappedAlphabet(),
            GappedAlphabet(),
            ReducedAlphabet("(AILMV)(NQST)(RHK)(DE)(FWY)CGP"),
        )
            for N = 1:3

                table = ContingencyTable(Float64, Val{N}, alphabet)

                fill!(table, AdditiveSmoothing(1.0))
                @test table[1] == 1.0
                @test getmarginals(table)[1] == length(alphabet)^(N - 1)
                @test sum(table) == length(alphabet)^N
                @test sum(getmarginals(table)) == float(N) * (length(alphabet)^N)
                @test gettotal(table) == float(length(alphabet))^N
            end
        end
    end

    @testset "Pseudocount" begin

        @test AdditiveSmoothing(1.0) == one(AdditiveSmoothing{Float64})
        @test AdditiveSmoothing(0.0) == zero(AdditiveSmoothing{Float64})

        for alphabet in (
            UngappedAlphabet(),
            GappedAlphabet(),
            ReducedAlphabet("(AILMV)(NQST)(RHK)(DE)(FWY)CGP"),
        )

            for N = 1:3

                table = ContingencyTable(Float64, Val{N}, alphabet)

                apply_pseudocount!(table, zero(AdditiveSmoothing{Float64}))

                @test sum(table.temporal) == 0.0
                @test sum(table) == 0.0
                @test sum(getmarginals(table)) == 0.0
                @test gettotal(table) == 0.0

                apply_pseudocount!(table, AdditiveSmoothing(1.0))

                @test table[1] == 1.0
                @test getmarginals(table)[1] == length(alphabet)^(N - 1)
                @test sum(table) == length(alphabet)^N
                @test sum(getmarginals(table)) == float(N) * (length(alphabet)^N)
                @test gettotal(table) == float(length(alphabet))^N

                table = ContingencyTable(Float64, Val{N}, alphabet)
                apply_pseudocount!(table, 1.0)

                @test table[1] == 1.0
                @test getmarginals(table)[1] == length(alphabet)^(N - 1)
                @test sum(table) == length(alphabet)^N
                @test sum(getmarginals(table)) == float(N) * (length(alphabet)^N)
                @test gettotal(table) == float(length(alphabet))^N
            end
        end
    end

end
