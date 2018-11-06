@testset "MSA Stats" begin

    msa_types = (
        Matrix{Residue},
        NamedResidueMatrix{Array{Residue,2}},
        MultipleSequenceAlignment,
        AnnotatedMultipleSequenceAlignment
        )

    pf09645_sto = joinpath(pwd(), "data", "PF09645_full.stockholm")
    gaoetal2011 = joinpath(pwd(), "data", "Gaoetal2011.fasta")

    gaoetal_msas = [ read(gaoetal2011, FASTA, T) for T in msa_types ]
    pfam_msas    = [ read(pf09645_sto, Stockholm, T) for T in msa_types ]

    @testset "Gaps" begin

        for msa in gaoetal_msas
            @test gapfraction(msa) == 0.0
            @test gapfraction(getsequence(msa,1)) == 0.0
            @test gapfraction(msa[1,:]) == 0.0

            @test residuefraction(msa) == 1.0

            @test sum(coverage(msa)) == 6.0
        end

        for msa in pfam_msas
            @test gapfraction(getsequence(msa,4)) == 0.0
            @test gapfraction(msa[:,1]) == 0.75
            @test gapfraction(msa[:,1]) == 0.75


            @test residuefraction(getsequence(msa, 4)) == 1.0
            @test residuefraction(msa[:,1]) == 0.25

            @test coverage(msa)[4] == 1.0
        end

    end
end
