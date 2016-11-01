@testset "Alphabet" begin

    @testset "Creation, Iteration and getindex" begin

        # Creation & Iteration
        @test [ i for i in UngappedAlphabet() ] == Int[ i for i in 1:20 ]
        @test [ i for i in GappedAlphabet() ]   == Int[ i for i in 1:21 ]
        @test [ i for i in ReducedAlphabet("(AILMV)(NQST)(RHK)(DE)(FWY)CGP") ]  ==
            Int[ i for i in 1:8 ]
    end

    @testset "in" begin

        # Creation & Iteration
        @test in(Residue('A'), UngappedAlphabet())
        @test in(Residue('A'), GappedAlphabet())
        @test in(Residue('A'), ReducedAlphabet("(AILMV)(NQST)(RHK)(DE)(FWY)CGP"))

        @test !in(GAP, UngappedAlphabet())
        @test  in(GAP, GappedAlphabet())
        @test !in(GAP, ReducedAlphabet("(AILMV)(NQST)(RHK)(DE)(FWY)CGP"))

        @test !in(XAA, UngappedAlphabet())
        @test !in(XAA, GappedAlphabet())
        @test !in(XAA, ReducedAlphabet("(AILMV)(NQST)(RHK)(DE)(FWY)CGP"))
    end

    @testset "Show" begin

        tmp = IOBuffer()
        show(tmp, UngappedAlphabet())
        @test takebuf_string(tmp) ==
        "MIToS.MSA.UngappedAlphabet of length 20. Residues : res\"ARNDCQEGHILKMFPSTWYV\""

        show(tmp, GappedAlphabet())
        @test takebuf_string(tmp) ==
        "MIToS.MSA.GappedAlphabet of length 21. Residues : res\"ARNDCQEGHILKMFPSTWYV-\""

        show(tmp, ReducedAlphabet("(AILMV)(NQST)(RHK)(DE)(FWY)CGP"))
        @test takebuf_string(tmp) ==
        "MIToS.MSA.ReducedAlphabet of length 8 : \"(AILMV)(NQST)(RHK)(DE)(FWY)CGP\""
    end
end
