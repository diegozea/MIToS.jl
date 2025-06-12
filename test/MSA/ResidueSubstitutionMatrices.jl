using MIToS.MSA.ResidueSubstitutionMatrices

@testset "ResidueSubstitutionMatrices" begin

    ab = GappedAlphabet()
    scores = Float64[i == j ? 1.0 : 0.0 for i = 1:length(ab), j = 1:length(ab)]

    m = ResidueSubstitutionMatrices.ResidueSubstitutionMatrix(scores)

    @test m[Residue('A'), Residue('A')] == 1.0
    @test m[Residue('A'), Residue('R')] == 0.0
    @test m[1, 1] == 1.0
    @test m[1, 2] == 0.0

    @test m[GAP, :] == scores[end, :]
    @test m[:, GAP] == scores[:, end]

    @test_throws ArgumentError ResidueSubstitutionMatrices.ResidueSubstitutionMatrix(
        ones(20, 21),
    )
    @test_throws ArgumentError ResidueSubstitutionMatrices.ResidueSubstitutionMatrix(
        ones(20, 20),
        ReducedAlphabet("AC"),
    )

    buf = IOBuffer()
    show(buf, MIME"text/plain"(), m)
    str = String(take!(buf))
    @test occursin("Named Matrix", str)
    @test occursin("A", str)
    @test !occursin("Alphabet:", str)

end
