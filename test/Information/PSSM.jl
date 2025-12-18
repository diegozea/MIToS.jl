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

    @testset "Alphabet controls gap counting" begin
        msa = Matrix{Residue}(undef, 3, 1)
        msa[:, 1] = Residue[Residue('A'), GAP, GAP]

        result = pssm(msa; alphabet = alphabet, background = uniform_background)

        a_index = alphabet[Residue('A')]
        c_index = alphabet[Residue('C')]

        @test isapprox(result.scores[a_index, 1], log2(1.0 / uniform_background[a_index]))
        @test result.scores[c_index, 1] == -Inf
    end

    @testset "GappedAlphabet counts gaps" begin
        gapped = GappedAlphabet()
        uniform_21 = fill(1 / length(gapped), length(gapped))
        msa = fill(GAP, 3, 1)

        result = pssm(msa; alphabet = gapped, background = uniform_21)

        gap_index = gapped[GAP]
        a_index = gapped[Residue('A')]
        q_gap = uniform_21[gap_index]

        @test isfinite(result.scores[gap_index, 1])
        @test isapprox(result.scores[gap_index, 1], log2(1.0 / q_gap))
        @test result.scores[a_index, 1] == -Inf
    end

    @testset "Zero probabilities" begin
        msa = Matrix{Residue}(undef, 2, 1)
        msa[:, 1] = Residue[Residue('A'), Residue('A')]

        c_index = alphabet[Residue('C')]
        result = pssm(msa; alphabet = alphabet, background = uniform_background)
        @test isinf(result.scores[c_index, 1]) && result.scores[c_index, 1] < 0
    end

    @testset "Zero background entries" begin
        q = copy(uniform_background)
        a_index = alphabet[Residue('A')]
        q[a_index] = 0.0
        q ./= sum(q)

        msa_a = fill(Residue('A'), 2, 1)
        result_inf = pssm(msa_a; alphabet = alphabet, background = q)
        @test isinf(result_inf.scores[a_index, 1]) && result_inf.scores[a_index, 1] > 0

        msa_c = fill(Residue('C'), 2, 1)
        result_nan = pssm(msa_c; alphabet = alphabet, background = q)
        @test isnan(result_nan.scores[a_index, 1])
    end

    @testset "Non-default base" begin
        msa = Matrix{Residue}(undef, 3, 1)
        msa[:, 1] = Residue[Residue('A'), Residue('A'), Residue('C')]

        result = pssm(msa; alphabet = alphabet, background = uniform_background, base = 10)

        p_a = 2 / 3
        q_a = uniform_background[alphabet[Residue('A')]]
        @test isapprox(result.scores[alphabet[Residue('A')], 1], log10(p_a / q_a))
    end

    @testset "All-gap column" begin
        msa = fill(GAP, 3, 1)
        result = pssm(msa; alphabet = alphabet, background = uniform_background)
        @test all(isnan, result.scores[:, 1])
    end

    @testset "dims validation" begin
        msa = reshape(res"AC", 1, :)
        @test_throws ArgumentError pssm(msa; dims = 1)
    end
end
