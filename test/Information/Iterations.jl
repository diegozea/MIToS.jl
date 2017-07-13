@testset "Iterations" begin

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
    @test isapprox(full(getarray(nmi)), result, rtol=1e-4)
end
