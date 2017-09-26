@testset "Information Measures" begin

    s = res"ARNDCQEGHILKMFPSTWYV-"
    r = reverse(s)
    g = res"GGGGGGGGGGGGGGGGGGGGG"
    Pg = probabilities(g)
    Pgg = probabilities(g,g)
    Pggg = probabilities(g,g,g)
    Ps = probabilities(s)
    Pss = probabilities(s,s)
    Psr = probabilities(s,r)
    Psss = probabilities(s,s,s)
    count_random = count(Residue[], Residue[], pseudocounts=AdditiveSmoothing(1.0))
    prob_random = probabilities(Residue[], Residue[], pseudocounts=AdditiveSmoothing(1.0))

    @testset "Entropy" begin

        @testset "H(X)" begin

            @test entropy(Pg) ≈ 0.0
            @test entropy(Ps) ≈ log(20.0)
            @test entropy(Ps, 2) ≈ log(2, 20.0)
        end

        @testset "Joint Entropy: H(X,Y)" begin

            @test entropy(Pgg) ≈ 0.0
            @test entropy(Pss) ≈ log(20.0)
            @test entropy(Psr) ≈ log(19.0)
            @test entropy(Pss, 2) ≈ log(2, 20.0)
            @test entropy(count_random) ≈ log(400.0)
            @test entropy(prob_random) ≈ log(400.0)
        end

        @testset "Joint Entropy: H(X,Y,Z)" begin

            @test entropy(Pggg) ≈ 0.0
            @test entropy(Psss) ≈ log(20)
        end

        @testset "Using counts" begin

            @test entropy(count(g,g)) ≈ entropy(Pgg)
            @test entropy(count(s,s)) ≈ entropy(Pss)
            @test entropy(count(s,r)) ≈ entropy(Psr)
        end
    end

    @testset "Kullback-Leibler" begin

        @test kullback_leibler(Ps, [1.0/20.0 for i in 1:20]) ≈ 0.0
        @test kullback_leibler(Ps, Ps) ≈ 0.0
        @test kullback_leibler(Ps, BLOSUM62_Pi) ≈ kullback_leibler(Ps)
        @test kullback_leibler(Ps, BLOSUM62_Pi) ≈ mapreduce(i->0.05*log(0.05/BLOSUM62_Pi[i]),+,1:20)
        @test kullback_leibler(Psr, Psr) ≈ 0.0
    end

    @testset "Mutual Information" begin

        @test mutual_information(Pgg) ≈ 0.0

        P = gettablearray(prob_random)
        Q = getmarginalsarray(prob_random)[:,1] * getmarginalsarray(prob_random)[:,2]' # 0.002500000000000001
        @test mutual_information(prob_random) ≈ sum(P .* log.(P./Q)) # So, it's != 0.0
        @test mutual_information(count_random) ≈ 0.0 # It's more accurate

        Pgs = probabilities(g,s)
        P = gettablearray(Pgs)
        Q = getmarginalsarray(Pgs)[:,1] * getmarginalsarray(Pgs)[:,2]' # 0.05000000000000002
        @test mutual_information(Pgs) ≈ mapreduce(x -> isnan(x) ? 0.0 : x, +, P.*log.(P./Q))
        @test mutual_information(count(g,s)) ≈ 0.0 # It's more accurate

        # MI(X,Y)=H(X)+H(Y)-H(X,Y)

        @test mutual_information(Psr) ≈ (marginal_entropy(Psr,1) +
                                         marginal_entropy(Psr,2) - entropy(Psr))

        @test mutual_information(Pss) ≈ (marginal_entropy(Pss,1) +
                                         marginal_entropy(Pss,2) - entropy(Pss))

        @test mutual_information(Psr,2) ≈ (marginal_entropy(Psr,1,2) +
                                           marginal_entropy(Psr,2,2) - entropy(Psr,2))

        @test mutual_information(Pss,20) ≈ (marginal_entropy(Pss,1,20) +
                                            marginal_entropy(Pss,2,20) - entropy(Pss,20))

        @test isapprox(mutual_information(prob_random), 0.0, atol=1e-15)
        @test mutual_information(count_random) ≈ 0.0
        @test isapprox(mutual_information(count_random), marginal_entropy(count_random,1) +
                       marginal_entropy(count_random,2) - entropy(count_random), atol=1e-13)

        @testset "Using counts" begin

            @test mutual_information(count_random) ≈ 0.0
            @test mutual_information(count(g,g)) ≈ mutual_information(Pgg)
            @test mutual_information(count(s,s)) ≈ mutual_information(Pss)
            @test mutual_information(count(s,r)) ≈ mutual_information(Psr)
        end

        @testset "MI(X,Y,Z)" begin

            @test mutual_information(Pggg) ≈ 0.0
            @test mutual_information(probabilities(g,s,r)) ≈ 0.0
            @test mutual_information(probabilities(s,s,r)) ≈ mutual_information(count(s,s,r))

            # MI(X,Y,Z) = H(X) + H(Y) + H(Z) - H(X,Y) - H(X,Z) - H(Y,Z) + H(X,Y,Z)
            @test mutual_information(Psss) ≈ mutual_information(count(s,s,s))
            @test mutual_information(Psss) ≈ (marginal_entropy(Psss, 1) +
                                              marginal_entropy(Psss, 2) +
                                              marginal_entropy(Psss, 3) -
                                              entropy(Pss) - entropy(Pss) -
                                              entropy(Pss) + entropy(Psss))

            # MI(X,Y,Z) <= min{ H(X,Y), H(X,Z), H(Y,Z) }
            @test mutual_information(Psss) ≈ mutual_information(Pss)
        end

        @testset "Pairwise Gap Percentage" begin

            @test gap_union_percentage(count(res"AA--", res"--AA", alphabet=GappedAlphabet())) ≈ 100.0
            @test gap_intersection_percentage(count(res"AA--", res"--AA", alphabet=GappedAlphabet())) ≈ 0.0

            @test gap_union_percentage(count(res"AA--", res"--AA", alphabet=GappedAlphabet(), weights=Weights([.25, .25, .25, .25]))) ≈ 100.0
            @test gap_intersection_percentage(count(res"AA--", res"--AA", alphabet=GappedAlphabet(), weights=Weights([.25, .25, .25, .25]))) ≈ 0.0

            @test gap_union_percentage(count(res"AAA-", res"AA--", alphabet=GappedAlphabet())) ≈ 50.0
            @test gap_intersection_percentage(count(res"AAA-", res"AA--", alphabet=GappedAlphabet())) ≈ 25.0

            @test gap_union_percentage(count(res"AAA-", res"AA--", alphabet=GappedAlphabet(), weights=Weights([.2, .2, .2, .4]))) ≈ 60.0
            @test gap_intersection_percentage(count(res"AAA-", res"AA--", alphabet=GappedAlphabet(), weights=Weights([.2, .2, .2, .4]))) ≈ 40.0
        end
    end
end
