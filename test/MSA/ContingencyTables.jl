@testset "ContingencyTables" begin

    @testset "Creation" begin

        for alphabet in (UngappedAlphabet(),
                         GappedAlphabet(),
                         ReducedAlphabet("(AILMV)(NQST)(RHK)(DE)(FWY)CGP"))
            for N in 1:3
                table = ContingencyTable(Float64, N, alphabet)

                @test size(table) == (Int[length(alphabet) for i in 1:N]...)
                @test length(table) == length(alphabet)^N
                @test sum(table) == 0.0
            end
        end
    end

    @testset "Indexing" begin

        for alphabet in (UngappedAlphabet(),
                         GappedAlphabet(),
                         ReducedAlphabet("(AILMV)(NQST)(RHK)(DE)(FWY)CGP"))

            table = ContingencyTable(Float64, 2, alphabet)

            if isa(get_alphabet(table), ReducedAlphabet)
                @test table[1,3] == 0.0
                table[1,3] = 10.0
                @test table[Residue('A'), Residue('R')] == 10.0
                table[Residue('A'), Residue('R')] = 20.0
                @test table["AILMV", "RHK"] == 20.0
                table["AILMV", "RHK"] = 30.0
                @test table[1,3] == 30.0
            else
                @test table[1,2] == 0.0
                table[1,2] = 10.0
                @test table[Residue('A'), Residue('R')] == 10.0
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

                table = ContingencyTable(Float64, N, alphabet)

                fill!(table.temporal, 1.0)
                @test sum(table.temporal) == 22.0^N
                @test sum(table) == 0.0
                @test sum(get_marginals(table)) == 0.0
                @test get_total(table) == 0.0

                MSA._update!(table)
                @test sum(table.temporal) == 22.0^N
                if isa(get_alphabet(table),ReducedAlphabet)
                    @test table[1] == 5.0^N
                    @test sum(table) == 20.0^N
                    @test sum(get_marginals(table)) == N*(20.0^N)
                    @test get_total(table) == 20.0^N
                    if N == 1
                        @test vec(array(get_marginals(table))) ==
                            vec(array(get_table(table)))
                    elseif N == 2
                        @test get_marginals(table)[1] ==
                            sum(5.0*n for n in [5,4,3,2,3,1,1,1])
                    elseif N == 2
                        @test get_marginals(table)[1] ==
                            sum(5.0*n*20.0 for n in [5,4,3,2,3,1,1,1])
                    end
                else
                    @test table[1] == 1.0
                    @test get_marginals(table)[1] == length(alphabet)^(N-1)
                    @test sum(table) == length(alphabet)^N
                    @test sum(get_marginals(table)) == float(N)*(length(alphabet)^N)
                    @test get_total(table) == float(length(alphabet))^N
                end
            end
        end
    end

    @testset "Fill & pseudocount" begin

        for alphabet in (UngappedAlphabet(),
                         GappedAlphabet(),
                         ReducedAlphabet("(AILMV)(NQST)(RHK)(DE)(FWY)CGP"))
            for N in 1:3

                table = ContingencyTable(Float64, N, alphabet)
                @test sum(table.temporal) == 0.0
                @test sum(table) == 0.0
                @test sum(get_marginals(table)) == 0.0
                @test get_total(table) == 0.0

                fill!(table, 1.0)
                @test table[1] == 1.0
                @test get_marginals(table)[1] == length(alphabet)^(N-1)
                @test sum(table) == length(alphabet)^N
                @test sum(get_marginals(table)) == float(N)*(length(alphabet)^N)
                @test get_total(table) == float(length(alphabet))^N
            end
        end
    end
end

