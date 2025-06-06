@testset "MSAQuality" begin
    import MIToS.MSA: _pairwise_score

    @testset "Sum Of Pairs" begin
        simple = read_file(joinpath(DATA, "simple.fasta"), FASTA)
        expected_simple =
            _pairwise_score(simple[1, :], simple[2, :]; gap_open = -10, gap_extend = -1)
        @test sum_of_pairs_score(simple) == expected_simple

        msa = permutedims(hcat(res"A-", res"AR", res"RA"))
        expected_default =
            _pairwise_score(msa[1, :], msa[2, :]; gap_open = -10, gap_extend = -1) +
            _pairwise_score(msa[1, :], msa[3, :]; gap_open = -10, gap_extend = -1) +
            _pairwise_score(msa[2, :], msa[3, :]; gap_open = -10, gap_extend = -1)
        @test sum_of_pairs_score(msa) == expected_default

        expected_custom =
            _pairwise_score(msa[1, :], msa[2, :]; gap_open = -5, gap_extend = -2) +
            _pairwise_score(msa[1, :], msa[3, :]; gap_open = -5, gap_extend = -2) +
            _pairwise_score(msa[2, :], msa[3, :]; gap_open = -5, gap_extend = -2)
        @test sum_of_pairs_score(msa; gap_open = -5, gap_extend = -2) == expected_custom

        custom_matrix = fill(-1, 22, 22)
        for i = 1:22
            custom_matrix[i, i] = 1
        end
        expected_matrix =
            _pairwise_score(
                msa[1, :],
                msa[2, :];
                gap_open = -5,
                gap_extend = -2,
                matrix = custom_matrix,
            ) +
            _pairwise_score(
                msa[1, :],
                msa[3, :];
                gap_open = -5,
                gap_extend = -2,
                matrix = custom_matrix,
            ) +
            _pairwise_score(
                msa[2, :],
                msa[3, :];
                gap_open = -5,
                gap_extend = -2,
                matrix = custom_matrix,
            )
        @test sum_of_pairs_score(
            msa;
            gap_open = -5,
            gap_extend = -2,
            matrix = custom_matrix,
        ) == expected_matrix

        g1 = res"-A-"
        g2 = res"A--"
        @test _pairwise_score(g1, g2; gap_open = -10, gap_extend = -1) == -20
    end
end
