@testset "PSSM" begin

    alphabet = UngappedAlphabet()
    uniform_background = fill(1 / 20, 20)

    @testset "Shape and ordering" begin
        msa = Matrix{Residue}(undef, 3, 2)
        msa[:, 1] = Residue[Residue('A'), Residue('C'), Residue('D')]
        msa[:, 2] = Residue[Residue('E'), Residue('F'), Residue('G')]

        result = pssm(msa; alphabet = alphabet, background = uniform_background)

        @test size(result.scores) == (length(alphabet), size(msa, 2))
        @test result.alphabet == alphabet
        @test result.background ≈ uniform_background
    end

    @testset "Log-odds correctness" begin
        msa = Matrix{Residue}(undef, 3, 1)
        msa[:, 1] = Residue[Residue('A'), Residue('A'), Residue('C')]

        result = pssm(msa; alphabet = alphabet, background = uniform_background, base = 2)

        p_a = 2 / 3
        p_c = 1 / 3
        q_a = uniform_background[alphabet[Residue('A')]]
        q_c = uniform_background[alphabet[Residue('C')]]

        @test isapprox(result.scores[alphabet[Residue('A')], 1], log2(p_a / q_a))
        @test isapprox(result.scores[alphabet[Residue('C')], 1], log2(p_c / q_c))
    end

    @testset "Zero policies" begin
        msa = Matrix{Residue}(undef, 2, 1)
        msa[:, 1] = Residue[Residue('A'), Residue('A')]

        neg_inf = pssm(msa; alphabet = alphabet, background = uniform_background, zero_policy = :negInf)
        clamp = pssm(msa; alphabet = alphabet, background = uniform_background, zero_policy = :clamp)

        c_index = alphabet[Residue('C')]
        @test isinf(neg_inf.scores[c_index, 1]) && neg_inf.scores[c_index, 1] < 0

        expected_clamp = log2(eps(Float64) / uniform_background[c_index])
        @test isapprox(clamp.scores[c_index, 1], expected_clamp)

        @test_throws DomainError pssm(
            msa;
            alphabet = alphabet,
            background = uniform_background,
            zero_policy = :error,
        )
    end

    @testset "All-gap column" begin
        msa = fill(GAP, 3, 1)
        result = pssm(
            msa;
            alphabet = alphabet,
            background = uniform_background,
            all_gap_value = 123.0,
        )

        @test all(result.scores[:, 1] .== 123.0)
    end

    @testset "dims validation" begin
        msa = reshape(res"AC", 1, :)
        @test_throws ArgumentError pssm(msa; dims = 1)
    end
end
