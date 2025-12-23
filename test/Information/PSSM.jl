@testset "PSSM" begin

    alphabet = UngappedAlphabet()
    uniform_background = fill(1 / 20, 20) # uniform distribution for the UngappedAlphabet

    @testset "Manual construction / empty matrices" begin
        pfm = PositionFrequencyMatrix(alphabet, 3)
        ppm = PositionSpecificProbabilityMatrix(alphabet, 3)
        pssm = PositionSpecificScoringMatrix(alphabet, 3, 2)

        @test Base.ismutabletype(typeof(pfm))
        @test Base.ismutabletype(typeof(ppm))
        @test Base.ismutabletype(typeof(pssm))

        @test all(==(0.0), pfm)
        @test all(isnan, ppm)
        @test all(==(0.0), pssm)
        @test pssm.base == 2

        pfm_zeros = zeros(PositionFrequencyMatrix, alphabet, 2)
        pssm_zeros = zeros(PositionSpecificScoringMatrix, alphabet, 2, 10)

        @test all(==(0.0), pfm_zeros)
        @test all(==(0.0), pssm_zeros)
        @test pssm_zeros.base == 10

        pssm[Residue('A'), 1] = 1.23
        @test pssm[Residue('A'), 1] == 1.23
        @test pssm["A", 1] == 1.23

        data = zeros(Float64, length(alphabet), 4)
        pssm2 = PositionSpecificScoringMatrix(data, alphabet, 2)
        pssm2[Residue('A'), 1] = 4.56
        @test pssm2[Residue('A'), 1] == 4.56
        @test gettablearray(pssm2)[alphabet[Residue('A')], 1] == 4.56

        pfm_named = PositionFrequencyMatrix(pfm.table, alphabet)
        @test pfm_named["A", 1] == pfm["A", 1]

        ppm_named = PositionSpecificProbabilityMatrix(ppm.table, alphabet)
        @test isnan(ppm_named["A", 1])

        @test_throws ArgumentError PositionFrequencyMatrix(pfm.table, GappedAlphabet())

        pfm_data = PositionFrequencyMatrix(zeros(Float64, length(alphabet), 2), alphabet)
        @test all(==(0.0), pfm_data)

        bad_data = zeros(Float64, length(alphabet) - 1, 2)
        @test_throws ErrorException PositionFrequencyMatrix(bad_data, alphabet)

        pssm_float = PositionSpecificScoringMatrix(alphabet, 2, 2.5)
        @test pssm_float.base == 2.5
        @test pssm_float.base isa Float64

        pssm_int = PositionSpecificScoringMatrix(alphabet, 2, 2)
        @test pssm_int.base === 2

        pssm_nat = PositionSpecificScoringMatrix(alphabet, 2, ℯ)
        @test pssm_nat.base === ℯ

        for base in (-1.0, 0.0, 1.0, NaN)
            @test_throws ArgumentError PositionSpecificScoringMatrix(alphabet, 2, base)
        end
    end

    @testset "Score sequence" begin
        pssm = PositionSpecificScoringMatrix(alphabet, 2)
        pssm["A", 1] = 1.0
        pssm["C", 2] = 2.0
        pssm["A", 2] = 0.5
        pssm["C", 1] = -1.0

        seq = Residue['A', 'C']
        pssm_score = score_sequence(pssm, seq)
        @test pssm_score isa ProfileScore
        @test pssm_score.score == 3.0
        @test pssm_score.used_positions == 2
        @test pssm_score.kind == :log_odds

        @test_throws ArgumentError score_sequence(pssm, Residue['A'])

        seqs_row = reshape(seq, 1, :)
        seqs_col = reshape(seq, :, 1)
        pssm_score_row = score_sequence(pssm, seqs_row)
        @test pssm_score_row isa ProfileScore
        @test pssm_score_row.score == 3.0
        @test pssm_score_row.used_positions == 2
        @test pssm_score_row.kind == :log_odds
        pssm_score_col = score_sequence(pssm, seqs_col)
        @test pssm_score_col isa ProfileScore
        @test pssm_score_col.score == 3.0
        @test pssm_score_col.used_positions == 2
        @test pssm_score_col.kind == :log_odds

        seqs_bad = fill(Residue('A'), 2, 3)
        @test_throws ArgumentError score_sequence(pssm, seqs_bad)

        seq_unknown = Residue['A', Residue('X')]
        score_unknown = @test_logs (:warn, r"Residue .* not in alphabet") score_sequence(
            pssm,
            seq_unknown,
        )
        @test score_unknown isa ProfileScore
        @test score_unknown.score == 1.0
        @test score_unknown.used_positions == 1
        @test score_unknown.kind == :log_odds

        seq_gap = Residue['A', GAP]
        score_gap = score_sequence(pssm, seq_gap)
        @test score_gap isa ProfileScore
        @test score_gap.score == 1.0
        @test score_gap.used_positions == 1
        @test score_gap.kind == :log_odds

        ppm = PositionSpecificProbabilityMatrix(alphabet, 2)
        ppm["A", 1] = 0.1
        ppm["C", 2] = 0.2
        ppm_score = score_sequence(ppm, seq)
        @test ppm_score isa ProfileScore
        @test isapprox(ppm_score.score, log(0.1) + log(0.2))
        @test ppm_score.base == ℯ
        @test isapprox(ppm_score.base^ppm_score.score, 0.02)
        @test ppm_score.used_positions == 2
        @test ppm_score.kind == :log_likelihood

        ppm_score_row = score_sequence(ppm, seqs_row)
        @test ppm_score_row isa ProfileScore
        @test isapprox(ppm_score_row.score, ppm_score.score)
        @test ppm_score_row.base == ppm_score.base
        @test ppm_score_row.used_positions == ppm_score.used_positions
        @test ppm_score_row.kind == :log_likelihood

        ppm_score_bits = score_sequence(ppm, seq; base = 2)
        @test ppm_score_bits isa ProfileScore
        @test isapprox(ppm_score_bits.score, ppm_score.score / log(2))
        @test ppm_score_bits.base == 2
        @test isapprox(ppm_score_bits.base^ppm_score_bits.score, 0.02)
        @test ppm_score_bits.kind == :log_likelihood

        seq_unknown_ppm = Residue['A', Residue('X')]
        ppm_unknown = @test_logs (:warn, r"Residue .* not in alphabet") score_sequence(
            ppm,
            seq_unknown_ppm,
        )
        @test ppm_unknown isa ProfileScore
        @test isapprox(ppm_unknown.score, log(0.1))
        @test ppm_unknown.base == ℯ
        @test isapprox(ppm_unknown.base^ppm_unknown.score, 0.1)
        @test ppm_unknown.used_positions == 1
        @test ppm_unknown.kind == :log_likelihood

        seq_gap_ppm = Residue['A', GAP]
        ppm_gap = score_sequence(ppm, seq_gap_ppm)
        @test ppm_gap isa ProfileScore
        @test isapprox(ppm_gap.score, log(0.1))
        @test ppm_gap.base == ℯ
        @test isapprox(ppm_gap.base^ppm_gap.score, 0.1)
        @test ppm_gap.used_positions == 1
        @test ppm_gap.kind == :log_likelihood

        ppm_zero = PositionSpecificProbabilityMatrix(alphabet, 2)
        ppm_zero["A", 1] = 0.1
        ppm_zero["C", 2] = 0.0
        zero_score = score_sequence(ppm_zero, seq)
        @test zero_score isa ProfileScore
        @test zero_score.score == -Inf
        @test zero_score.base == ℯ
        @test zero_score.used_positions == 2
        @test zero_score.kind == :log_likelihood

        @testset "MSA and sequence types" begin
            seq_vec = Residue['A', 'C', 'D']
            msa_matrix = permutedims(seq_vec)

            msa = MultipleSequenceAlignment(msa_matrix)
            annot_msa = AnnotatedMultipleSequenceAlignment(msa_matrix)

            pssm_types = PositionSpecificScoringMatrix(alphabet, 3)
            pssm_types["A", 1] = 0.5
            pssm_types["C", 2] = 1.0
            pssm_types["D", 3] = -0.25

            expected = score_sequence(pssm_types, seq_vec)
            @test expected.used_positions == 3
            @test expected.kind == :log_odds

            aligned_seq = getsequence(msa, 1)
            annot_aligned_seq = getsequence(annot_msa, 1)
            unaligned_seq = AnnotatedSequence(annot_aligned_seq)

            @test score_sequence(pssm_types, msa) == expected
            @test score_sequence(pssm_types, annot_msa) == expected
            @test score_sequence(pssm_types, aligned_seq) == expected
            @test score_sequence(pssm_types, annot_aligned_seq) == expected
            @test score_sequence(pssm_types, unaligned_seq) == expected
        end
    end

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
