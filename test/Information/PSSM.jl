@testset "PSSM" begin

    alphabet = UngappedAlphabet()
    uniform_background = fill(1 / 20, 20) # uniform distribution for the UngappedAlphabet

    @testset "Shape, ordering and values" begin
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
        @test length(result) == length(alphabet) * size(msa, 2)
        @test result.alphabet == alphabet
        @test eltype(result) === Float64
        @test sum(isfinite.(result), dims = 1) == [3 3] # 3 different residues per column
    end

    @testset "Frequencies" begin
        msa = Residue[
            'A' 'C'
            'A' GAP
            'C' 'C'
        ]

        pfm = position_frequency_matrix(msa; alphabet = alphabet)

        a_index = alphabet[Residue('A')] # 1
        c_index = alphabet[Residue('C')] # 5

        @test pfm.table[a_index, 1] == 2
        @test pfm.table[c_index, 1] == 1
        @test pfm.table[c_index, 2] == 2
        @test pfm.table[a_index, 2] == 0
        @test sum(pfm, dims = 1) == [3 2] # gaps not counted

        gapped = GappedAlphabet()
        gpfm = position_frequency_matrix(msa; alphabet = gapped)
        gap_index = gapped[GAP] # 21

        @test gpfm.table[gap_index, 1] == 0
        @test gpfm.table[gap_index, 2] == 1
        @test sum(gpfm, dims = 1) == [3 3] # gaps counted
    end

    @testset "Probabilities" begin
        msa = Residue[
            'A' 'C'
            'A' GAP
            'C' 'C'
        ]

        ppm = position_specific_probability_matrix(msa; alphabet = alphabet)

        # P should sum to 1 per column
        @test isapprox(sum(ppm.table[:, 1]), 1.0)
        @test isapprox(sum(ppm.table[:, 2]), 1.0)

        msa_gap = fill(GAP, 3, 1) # all-gap column
        ppm_gap = position_specific_probability_matrix(msa_gap; alphabet = alphabet)
        @test all(isnan, ppm_gap.table[:, 1])

        # same result when computed from frequencies than directly from an MSA
        ppm_from_freqs = position_specific_probability_matrix(
            position_frequency_matrix(msa; alphabet = alphabet),
        )
        @test isapprox(gettablearray(ppm), gettablearray(ppm_from_freqs))
    end

    @testset "Scoring overloads and non-mutating semantics" begin
        msa = reshape(Residue['A', 'A', 'C'], 3, 1)

        pssm_msa = position_specific_scoring_matrix(
            msa;
            alphabet = alphabet,
            background = uniform_background,
            base = 2,
        )
        pssm_freqs = position_specific_scoring_matrix(
            position_frequency_matrix(msa; alphabet = alphabet);
            background = uniform_background,
            base = 2,
        )
        pssm_probs = position_specific_scoring_matrix(
            position_specific_probability_matrix(msa; alphabet = alphabet);
            background = uniform_background,
            base = 2,
        )

        @test isapprox(gettablearray(pssm_msa), gettablearray(pssm_freqs))
        @test isapprox(gettablearray(pssm_msa), gettablearray(pssm_probs))

        @testset "non-mutating behavior" begin
            pfm_original = position_frequency_matrix(msa; alphabet = alphabet)
            pfm_backup = deepcopy(pfm_original)
            position_specific_probability_matrix(pfm_original)
            @test gettablearray(pfm_original) == gettablearray(pfm_backup)

            ppm_original = position_specific_probability_matrix(msa; alphabet = alphabet)
            ppm_backup = deepcopy(ppm_original)
            position_specific_scoring_matrix(
                ppm_original;
                background = uniform_background,
                base = 2,
            )
            @test gettablearray(ppm_original) == gettablearray(ppm_backup)
        end
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
        q_a = 0.05 # uniform background for the UngappedAlphabet
        q_c = 0.05

        @test isapprox(result.table[alphabet[Residue('A')], 1], log2(p_a / q_a))
        @test isapprox(result.table[alphabet[Residue('C')], 1], log2(p_c / q_c))
        @test isinf(result.table["R", 1]) # arginine not present in column
    end

    @testset "Alphabet controls gap counting" begin
        msa = reshape(Residue['A', GAP, GAP], 3, 1)

        result = position_specific_scoring_matrix(
            msa;
            alphabet = alphabet,
            background = uniform_background,
        )

        # UngappedAlphabet, so gaps not counted
        @test isapprox(result["A", 1], log(1.0 / 0.05)) # uniform background
        @test result["C", 1] == -Inf
        @test result[1, 1] == result["A", 1] # index access works the same
        @test result[1, 1] == result[Residue('A'), 1] # Residue access works the same
        @test_throws BoundsError result[21, 1] # there is no value for gaps
        @test_throws BoundsError result[GAP, 1] # there is no value for gaps
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
