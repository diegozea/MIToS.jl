@testset "Alphabet" begin

    @testset "Creation, Iteration and getindex" begin

        # Creation & Iteration
        @test [ i for i in UngappedAlphabet() ] == Int[ i for i in 1:20 ]
        @test [ i for i in GappedAlphabet() ]   == Int[ i for i in 1:21 ]
        @test [ i for i in ReducedAlphabet("(AILMV)(NQST)(RHK)(DE)(FWY)CGP") ]  ==
           Int[ i for i in 1:8 ]
    end

    @testset "getindex and names" begin

        reduced = ReducedAlphabet("(AILMV)(NQST)(RHK)(DE)(FWY)CGP")
        strings = ["AILMV", "NQST", "RHK", "DE", "FWY", "C", "G", "P"]

        @test length(reduced) == length(strings)
        @test names(reduced) == strings

        for i in 1:length(reduced)
            @test reduced[strings[i]] == i
        end

        strings_gapped = ["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S",
            "T","W","Y","V","-"]

        @test names(GappedAlphabet()) == strings_gapped
        @test names(UngappedAlphabet()) == strings_gapped[1:end-1]

        for i in 1:20
            @test UngappedAlphabet()[strings_gapped[i]] == i
            @test GappedAlphabet()[strings_gapped[i]] == i
        end

        @test GappedAlphabet()["-"] == 21
        @test UngappedAlphabet()["-"] == 22
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
        @test String(take!(tmp)) ==
        "UngappedAlphabet of length 20. Residues : res\"ARNDCQEGHILKMFPSTWYV\""

        show(tmp, GappedAlphabet())
        @test String(take!(tmp)) ==
        "GappedAlphabet of length 21. Residues : res\"ARNDCQEGHILKMFPSTWYV-\""

        show(tmp, ReducedAlphabet("(AILMV)(NQST)(RHK)(DE)(FWY)CGP"))
        @test String(take!(tmp)) ==
        "ReducedAlphabet of length 8 : \"(AILMV)(NQST)(RHK)(DE)(FWY)CGP\""
    end
end
