@testset "MSAQuality" begin
    import MIToS.MSA: _pairwise_score

    @testset "Sum Of Pairs" begin
        # simple.fasta : 
        # AR
        # RA
        simple = read_file(joinpath(DATA, "simple.fasta"), FASTA) 
        expected_simple = MSA.BLOSUM62[Int(Residue('A')), Int(Residue('R'))] +
                          MSA.BLOSUM62[Int(Residue('R')), Int(Residue('A'))]
        @test _pairwise_score(simple[1, :], simple[2, :]; gap_open = -10, gap_extend = -1) == expected_simple
        @test sum_of_pairs_score(simple) == expected_simple

        msa = Residue[
            'A'  '-'  '-'  'A'
            'A'  'R'  'A'  'V'
            '-'  'R'  'A'  '-'
        ]
        
        # BLOSUM62 values
        # A A =  4
        # A R = -1
        # R R =  5
        # A V =  0
        for (go, ge) in [(-10, -1), (-5, -2)] # (-10, -1) is the default
            # 1-2, 1-3, 2-3
            expected = (4 + go + ge + 0) + 
                       (go + go + ge + go) + 
                       (go + 5 + 4 + go)
            observed = _pairwise_score(msa[1, :], msa[2, :]; gap_open = -go, gap_extend = ge) +
                _pairwise_score(msa[1, :], msa[3, :]; gap_open = -go, gap_extend = ge) +
                _pairwise_score(msa[2, :], msa[3, :]; gap_open = -go, gap_extend = ge)
            @test observed == expected
            @test sum_of_pairs_score(msa; gap_open = go, gap_extend = ge) == expected
        end
        # Test with multiple gaps
        g1 = res"-A-"
        g2 = res"A--"
        @test _pairwise_score(g1, g2) == -20 # gap_open=-10
    end
end
