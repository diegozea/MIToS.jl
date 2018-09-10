@testset "ContingencyTables" begin

    @testset "Creation and Getters" begin

        for alphabet in (UngappedAlphabet(),
                         GappedAlphabet(),
                         ReducedAlphabet("(AILMV)(NQST)(RHK)(DE)(FWY)CGP"))
            for N in 1:3
                table = ContingencyTable(Float64, Val{N}, alphabet) # zeros in MIToS 1.0

                @test size(table) == (Int[length(alphabet) for i in 1:N]...,)
                @test length(table) == length(alphabet)^N
                @test length(getmarginals(table)) == length(alphabet)* N
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
            for alphabet in (UngappedAlphabet(),
                            GappedAlphabet(),
                            ReducedAlphabet("(AILMV)(NQST)(RHK)(DE)(FWY)CGP"))
                for N in 1:3
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

    @testset "Indexing" begin

        for alphabet in (UngappedAlphabet(),
                         GappedAlphabet(),
                         ReducedAlphabet("(AILMV)(NQST)(RHK)(DE)(FWY)CGP"))

            table = ContingencyTable(Float64, Val{2}, alphabet)

            if isa(getalphabet(table), ReducedAlphabet)
                @test table[1,3] == 0.0
                table[1,3] = 100.0
                i = 2 * length(getalphabet(table)) + 1 # using one index
                @test table[i] == 100.0
                table[i] = 10.0
                @test table[Residue('A'), Residue('R')] == 10.0 # using two indices
                table[Residue('A'), Residue('R')] = 20.0
                @test table["AILMV", "RHK"] == 20.0
                table["AILMV", "RHK"] = 30.0
                @test table[1,3] == 30.0
            else
                @test table[1,2] == 0.0
                table[1,2] = 100.0
                i = length(getalphabet(table)) + 1 # using one index
                @test table[i] == 100.0
                table[i] = 10.0
                @test table[Residue('A'), Residue('R')] == 10.0 # using two indices
                table[Residue('A'), Residue('R')] = 20.0
                @test table["A", "R"] == 20.0
                table["A", "R"] = 30.0
                @test table[1,2] == 30.0
             end
        end
    end

    @testset "Update" begin

        for alphabet in (UngappedAlphabet(),
                         GappedAlphabet(),
                         ReducedAlphabet("(AILMV)(NQST)(RHK)(DE)(FWY)CGP"))
            for N in 1:3

                table = ContingencyTable(Float64, Val{N}, alphabet)

                fill!(table.temporal, 1.0)
                @test sum(table.temporal) == 22.0^N
                @test sum(table) == 0.0
                @test sum(getmarginals(table)) == 0.0
                @test gettotal(table) == 0.0

                Information._update!(table)
                @test sum(table.temporal) == 22.0^N
                if isa(getalphabet(table),ReducedAlphabet)
                    @test table[1] == 5.0^N
                    @test sum(table) == 20.0^N
                    @test sum(getmarginals(table)) == N*(20.0^N)
                    @test gettotal(table) == 20.0^N
                    if N == 1
                        @test vec(getarray(getmarginals(table))) ==
                            vec(getarray(gettable(table)))
                    elseif N == 2
                        @test getmarginals(table)[1] ==
                            sum(5.0*n for n in [5,4,3,2,3,1,1,1])
                    elseif N == 2
                        @test getmarginals(table)[1] ==
                            sum(5.0*n*20.0 for n in [5,4,3,2,3,1,1,1])
                    end
                else
                    @test table[1] == 1.0
                    @test getmarginals(table)[1] == length(alphabet)^(N-1)
                    @test sum(table) == length(alphabet)^N
                    @test sum(getmarginals(table)) == float(N)*(length(alphabet)^N)
                    @test gettotal(table) == float(length(alphabet))^N
                end
            end
        end
    end

    @testset "Fill" begin

        for alphabet in (UngappedAlphabet(),
                         GappedAlphabet(),
                         ReducedAlphabet("(AILMV)(NQST)(RHK)(DE)(FWY)CGP"))
            for N in 1:3

                table = ContingencyTable(Float64, Val{N}, alphabet)

                fill!(table, AdditiveSmoothing(1.0))
                @test table[1] == 1.0
                @test getmarginals(table)[1] == length(alphabet)^(N-1)
                @test sum(table) == length(alphabet)^N
                @test sum(getmarginals(table)) == float(N)*(length(alphabet)^N)
                @test gettotal(table) == float(length(alphabet))^N
            end
        end
    end

    @testset "Pseudocount" begin

        @test AdditiveSmoothing(1.0) ==  one(AdditiveSmoothing{Float64})
        @test AdditiveSmoothing(0.0) == zero(AdditiveSmoothing{Float64})

        for alphabet in (UngappedAlphabet(),
                         GappedAlphabet(),
                         ReducedAlphabet("(AILMV)(NQST)(RHK)(DE)(FWY)CGP"))

            for N in 1:3

                table = ContingencyTable(Float64, Val{N}, alphabet)

                apply_pseudocount!(table, zero(AdditiveSmoothing{Float64}))

                @test sum(table.temporal) == 0.0
                @test sum(table) == 0.0
                @test sum(getmarginals(table)) == 0.0
                @test gettotal(table) == 0.0

                apply_pseudocount!(table, AdditiveSmoothing(1.0))

                @test table[1] == 1.0
                @test getmarginals(table)[1] == length(alphabet)^(N-1)
                @test sum(table) == length(alphabet)^N
                @test sum(getmarginals(table)) == float(N)*(length(alphabet)^N)
                @test gettotal(table) == float(length(alphabet))^N

                table = ContingencyTable(Float64, Val{N}, alphabet)
                apply_pseudocount!(table, 1.0)

                @test table[1] == 1.0
                @test getmarginals(table)[1] == length(alphabet)^(N-1)
                @test sum(table) == length(alphabet)^N
                @test sum(getmarginals(table)) == float(N)*(length(alphabet)^N)
                @test gettotal(table) == float(length(alphabet))^N
            end
        end
    end

end
