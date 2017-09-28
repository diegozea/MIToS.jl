@testset "Iterations" begin

    @testset "NMI" begin
        # This is the example of MI(X, Y)/H(X, Y) from:
        #
        # Gao, H., Dou, Y., Yang, J., & Wang, J. (2011).
        # New methods to measure residues coevolution in proteins.
        # BMC bioinformatics, 12(1), 206.

        aln = read(joinpath(pwd(), "data", "Gaoetal2011.fasta"), FASTA)
        result = Float64[ 0     0     0     0     0     0
                          0     0     0     0     0     0
                          0     0     0     1     1     0.296
                          0     0     1     0     1     0.296
                          0     0     1     1     0     0.296
                          0     0     0.296 0.296 0.296 0 ]

        nmi = mapcolpairfreq!(normalized_mutual_information, aln,
                              Counts(ContingencyTable(Float64,Val{2},UngappedAlphabet())),
                              Val{false})
        nmi_mat = convert(Matrix{Float64}, getarray(nmi))
        @test isapprox(nmi_mat, result, rtol=1e-4)

        nmi_t = mapseqpairfreq!(normalized_mutual_information, transpose(aln),
                                Counts(ContingencyTable(Float64,Val{2},UngappedAlphabet())),
                                Val{false})
        @test nmi_mat == convert(Matrix{Float64},getarray(nmi_t))
    end

    @testset "Gaps" begin

        function _gaps(table::Union{Probabilities{Float64,1,GappedAlphabet},
                                    Counts{Float64,1,GappedAlphabet}})
            table[GAP]
        end

        table = ContingencyTable(Float64, Val{1}, GappedAlphabet())

        gaps = read(joinpath(pwd(), "data", "gaps.txt"), Raw)

        # THAYQAIHQV 0
        # THAYQAIHQ- 0.1
        # THAYQAIH-- 0.2
        # THAYQAI--- 0.3
        # THAYQA---- 0.4
        # THAYQ----- 0.5
        # THAY------ 0.6
        # THA------- 0.7
        # TH-------- 0.8
        # T--------- 0.9

        ngaps = Float64[0,1,2,3,4,5,6,7,8,9]

        colcount = mapcolfreq!(_gaps, gaps, Counts(table))
        @test all((vec(getarray(colcount)) .- ngaps) .== 0.0)

        colfract = mapcolfreq!(_gaps, gaps, Probabilities(table))
        @test all((vec(getarray(colfract)) .- ngaps./10.0) .== 0.0)

        seqcount = mapseqfreq!(_gaps, gaps, Counts(table))
        @test all((vec(getarray(seqcount)) .- ngaps) .== 0.0)

        seqfract = mapseqfreq!(_gaps, gaps, Probabilities(table))
        @test all((vec(getarray(seqfract)) .- ngaps./10.0) .== 0.0)
    end
end
