using Pkg
Pkg.add("BioAlignments")
using MIToS.MSA.ResidueSubstitutionMatrices
using MIToS.MSA
using BioAlignments

@testset "BioAlignmentsExt" begin
    mito = MIToS.MSA.ResidueSubstitutionMatrices.BLOSUM62
    bio = BioAlignments.BLOSUM62
    conv = convert(ResidueSubstitutionMatrix, bio)
    @test conv.scores == mito.scores
    @test typeof(conv.alphabet) == typeof(mito.alphabet)
end
