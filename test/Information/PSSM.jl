@testset "PSSM" begin

    alphabet = UngappedAlphabet()
    uniform_background = fill(1 / 20, 20)

    @testset "Shape and ordering" begin
        msa = Residue[
            'A' 'E'
            'C' 'F'
            'D' 'G'
        ]

        result = position_specific_scoring_matrix(
            msa;
            alphabet = alphabet,
            background = uniform_background,
        )

        @test size(result.table) == (length(alphabet), size(msa, 2))
        @test result.alphabet == alphabet
    end

    @testset "Log-odds correctness" begin
        msa = reshape(Residue['A', 'A', 'C'], 3, 1)

        result = position_specific_scoring_matrix(
            msa;
            alphabet = alphabet,
            background = uniform_background,
            base = 2,
        )

        p_a = 2 / 3
        p_c = 1 / 3
        q_a = uniform_background[alphabet[Residue('A')]]
        q_c = uniform_background[alphabet[Residue('C')]]

        @test isapprox(result.table[alphabet[Residue('A')], 1], log2(p_a / q_a))
        @test isapprox(result.table[alphabet[Residue('C')], 1], log2(p_c / q_c))
    end

    @testset "Alphabet controls gap counting" begin
        msa = reshape(Residue['A', GAP, GAP], 3, 1)

        result = position_specific_scoring_matrix(
            msa;
            alphabet = alphabet,
            background = uniform_background,
        )

        a_index = alphabet[Residue('A')]
        c_index = alphabet[Residue('C')]

        @test isapprox(result.table[a_index, 1], log(1.0 / uniform_background[a_index]))
        @test result.table[c_index, 1] == -Inf
    end

    @testset "GappedAlphabet counts gaps" begin
        gapped = GappedAlphabet()
        uniform_21 = fill(1 / length(gapped), length(gapped))
        msa = fill(GAP, 3, 1)

        result = position_specific_scoring_matrix(
            msa;
            alphabet = gapped,
            background = uniform_21,
        )

        gap_index = gapped[GAP]
        a_index = gapped[Residue('A')]
        q_gap = uniform_21[gap_index]

        @test isfinite(result.table[gap_index, 1])
        @test isapprox(result.table[gap_index, 1], log(1.0 / q_gap))
        @test result.table[a_index, 1] == -Inf
    end

    @testset "Zero probabilities" begin
        msa = reshape(Residue['A', 'A'], 2, 1)

        c_index = alphabet[Residue('C')]
        result = position_specific_scoring_matrix(
            msa;
            alphabet = alphabet,
            background = uniform_background,
        )
        @test isinf(result.table[c_index, 1]) && result.table[c_index, 1] < 0
    end

    @testset "Zero background entries" begin
        q = copy(uniform_background)
        a_index = alphabet[Residue('A')]
        q[a_index] = 0.0
        q ./= sum(q)

        msa_a = fill(Residue('A'), 2, 1)
        result_inf =
            position_specific_scoring_matrix(msa_a; alphabet = alphabet, background = q)
        @test isinf(result_inf.table[a_index, 1]) && result_inf.table[a_index, 1] > 0

        msa_c = fill(Residue('C'), 2, 1)
        result_nan =
            position_specific_scoring_matrix(msa_c; alphabet = alphabet, background = q)
        @test isnan(result_nan.table[a_index, 1])
    end

    @testset "Non-default base" begin
        msa = reshape(Residue['A', 'A', 'C'], 3, 1)

        result = position_specific_scoring_matrix(
            msa;
            alphabet = alphabet,
            background = uniform_background,
            base = 10,
        )

        p_a = 2 / 3
        q_a = uniform_background[alphabet[Residue('A')]]
        @test isapprox(result.table[alphabet[Residue('A')], 1], log10(p_a / q_a))
    end

    @testset "Invalid base" begin
        msa = reshape(Residue['A'], 1, 1)

        for base in [-1.0, 0.0, 1.0]
            @test_throws ArgumentError position_specific_scoring_matrix(
                msa;
                alphabet = alphabet,
                background = uniform_background,
                base = base,
            )
        end
    end

    @testset "All-gap column" begin
        msa = fill(GAP, 3, 1)
        result = position_specific_scoring_matrix(
            msa;
            alphabet = alphabet,
            background = uniform_background,
        )
        @test all(isnan, result.table[:, 1])
    end

    @testset "Background validation" begin
        msa = permutedims(hcat(res"AC", res"AD")) # 2 sequences, 2 columns

        @testset "Length mismatch" begin
            background = fill(1 / 19, 19)
            @test_throws ArgumentError position_specific_scoring_matrix(
                msa;
                alphabet = alphabet,
                background = background,
            )
        end

        @testset "Zero sum" begin
            background = zeros(20)
            @test_throws DomainError position_specific_scoring_matrix(
                msa;
                alphabet = alphabet,
                background = background,
            )
        end

        @testset "Non-finite" begin
            background = fill(1 / 20, 20)
            background[1] = NaN
            @test_throws DomainError position_specific_scoring_matrix(
                msa;
                alphabet = alphabet,
                background = background,
            )
        end

        @testset "Negative" begin
            background = fill(1 / 20, 20)
            background[1] = -0.1
            @test_throws DomainError position_specific_scoring_matrix(
                msa;
                alphabet = alphabet,
                background = background,
            )
        end
    end
end
