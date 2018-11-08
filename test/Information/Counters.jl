@testset "Counters" begin

    @testset "count and probabilities" begin
        seq = res"ARNDCQEGHILKMFPSTWYV-"

        for alphabet in (UngappedAlphabet(),
                         GappedAlphabet(),
                         ReducedAlphabet("(AILMV)(NQST)(RHK)(DE)(FWY)CGP"))
            for N in 1:3
                seqs = ((seq for i in 1:N)...,)::NTuple{N,Vector{Residue}}

                table = ContingencyTable(Float64, Val{N}, alphabet) # zeros in MIToS 1.0
                if N == 1
                    @test table[1] == 0.0
                elseif N == 2
                    @test table[1,1] == 0.0
                else
                    @test table[1,1,1] == 0.0
                end
                @test getmarginals(table)[1,1] == 0.0
                @test gettotal(table) == 0.0

                count!(table, NoClustering(), NoPseudocount(), seqs...)

                @test table == count(seqs..., alphabet=alphabet)

                if isa(alphabet, ReducedAlphabet)
                    @test table[1] == 5.0
                    @test getmarginals(table)[1] == 5.0
                    @test gettotal(table) ≈ 20.0 # Reduced alphabet without gap
                else
                    @test table[1] == 1.0
                    @test getmarginals(table)[1] == 1.0
                    @test gettotal(table) ≈ length(alphabet)
                    if N == 2
                        len = length(alphabet)
                        @test gettablearray(table) == Matrix{Float64}(I, len, len)
                        @test getmarginalsarray(table)[:,1] == [1.0 for i in 1:len]
                    end
                end

                normalize!(table)
                @test table == probabilities(seqs..., alphabet=alphabet)
            end
        end

        @testset "Using clustering" begin

            clusters = Clusters([21], ones(Int,21), Weights(ones(21)/21))

            for alphabet in (UngappedAlphabet(),
                            GappedAlphabet(),
                            ReducedAlphabet("(AILMV)(NQST)(RHK)(DE)(FWY)CGP"))

                for N in 1:3
                    seqs = ((seq for i in 1:N)...,)::NTuple{N,Vector{Residue}}

                    table = ContingencyTable(Float64, Val{N}, alphabet)

                    count!(table, clusters, NoPseudocount(), seqs...)

                    @test table == count(seqs..., alphabet=alphabet, weights=clusters)

                    len = length(alphabet)
                    if isa(alphabet, ReducedAlphabet)
                        @test table[1] == 5.0/21
                        @test getmarginals(table)[1] == 5.0/21
                        @test gettotal(table) ≈ (1.0/21) * 20.0 # ReducedAlphabet without gap
                    else
                        @test table[1] == 1.0/21
                        @test getmarginals(table)[1] == 1.0/21
                        @test gettotal(table) ≈ (1.0/21) * len
                        if N == 2
                            @test gettablearray(table) == (1.0/21) .* Matrix{Float64}(I, len, len)
                            @test getmarginalsarray(table)[:,1] == [1.0/21 for i in 1:len]
                        end
                    end

                    normalize!(table)
                    @test table == probabilities(seqs..., alphabet=alphabet, weights=clusters)
                end
            end
        end

        @testset "Using pseudocount" begin

            for alphabet in (UngappedAlphabet(),
                            GappedAlphabet(),
                            ReducedAlphabet("(AILMV)(NQST)(RHK)(DE)(FWY)CGP"))

                for N in 1:3
                    seqs = ((seq for i in 1:N)...,)::NTuple{N,Vector{Residue}}

                    table = ContingencyTable(Float64, Val{N}, alphabet)

                    count!(table, NoClustering(), AdditiveSmoothing(1.0), seqs...)

                    @test table == count(seqs..., alphabet = alphabet,
                                                  pseudocounts = AdditiveSmoothing(1.0))

                    len = Float64(length(alphabet))
                    if isa(alphabet, ReducedAlphabet)
                        @test table[1] == 5.0 + 1.0
                        @test getmarginals(table)[1] == 5. + (N==1 ? 1. : N==2 ? len : len^2)
                        @test gettotal(table) ≈ length(gettable(table)) + 20.0 # without gap
                    else
                        @test table[1] == 2.0
                        @test getmarginals(table)[1] == 1. + (N==1 ? 1. : N==2 ? len : len^2)
                        @test gettotal(table) ≈ len + length(gettable(table))
                        if N == 2
                            @test gettablearray(table) == Matrix{Float64}(I, Int(len), Int(len)) .+ 1.0
                            @test getmarginalsarray(table)[:,1] == [ len + 1.0 for i in 1:len ]
                        end
                    end

                    normalize!(table)
                    @test table == probabilities(seqs..., alphabet = alphabet,
                                                 pseudocounts = AdditiveSmoothing(1.0))
                end
            end
        end

        @testset "pseudofrequencies" begin

            table = ContingencyTable(Float64, Val{2}, UngappedAlphabet())

            probabilities!(table, NoClustering(), NoPseudocount(),
                                  BLOSUM_Pseudofrequencies(1.,0.), seq, seq)

            @test table == probabilities(seq, seq,
                                         pseudofrequencies = BLOSUM_Pseudofrequencies(1.,0.))
            @test table == probabilities(seq, seq)

            @test table[1] == 1.0/20.
            @test getmarginals(table)[1] == 1.0/20.
            @test gettotal(table) ≈ 1.0
            @test gettablearray(table) == Matrix{Float64}(I, 20, 20) ./ 20.0

            @test table != probabilities(seq, seq,
                                         pseudofrequencies = BLOSUM_Pseudofrequencies(0.,1.))
        end

        @testset "BigFloat" begin
            table = ContingencyTable(BigFloat, Val{2}, UngappedAlphabet())
            clusters = Clusters([21], ones(Int,21), Weights(ones(21)/21))
            probabilities!(table, clusters, AdditiveSmoothing(one(BigFloat)),
                           BLOSUM_Pseudofrequencies(1.0,1.0), seq, seq)

            @test sum(table) ≈ one(BigFloat)
            @test eltype(table) == BigFloat
            @test isa(table[1], BigFloat)
        end

        @testset "delete_dimensions" begin
            Pxyz = probabilities(seq, seq, seq)
            @test delete_dimensions(Pxyz, 3) == probabilities(seq, seq)
            @test delete_dimensions(Pxyz, 3, 2) == probabilities(seq)
            Pxy = delete_dimensions(Pxyz, 3)
            @test delete_dimensions!(Pxy, Pxyz, 1) == probabilities(seq, seq)
            @test sum(Pxy) ≈ 1.0
            Nxyz = count(seq, seq, seq)
            Nxy = delete_dimensions(Nxyz, 3)
            @test Nxy == count(seq, seq)
            @test delete_dimensions(Nxyz, 3, 2) == count(seq)
            @test delete_dimensions!(Nxy, Nxyz, 1) == count(seq, seq)
        end
    end
end
