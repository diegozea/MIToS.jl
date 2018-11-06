@testset "Shuffle" begin

    msa_types = (
        Matrix{Residue},
        NamedResidueMatrix{Array{Residue,2}},
        MultipleSequenceAlignment,
        AnnotatedMultipleSequenceAlignment
        )

    N = length(msa_types)

    pf09645_sto = joinpath(pwd(), "data", "PF09645_full.stockholm")

    msas = [ read(pf09645_sto, Stockholm, T) for T in msa_types ]
    gaps = [ msa .== GAP for msa in msas ]
    lcol = [ mean(msa .== Residue('L'), dims=1) for msa in msas ]
    lseq = [ mean(msa .== Residue('L'), dims=2) for msa in msas ]

    Random.seed!(42)

    @testset "General" begin

        msa = msas[1]

        for dim in [1,2]

            aln = copy(msa)
            shuffle!(aln, dim)
            @test aln != msa
            @test (aln .== GAP) == gaps[1] # default: fixed gaps

            aln = copy(msa)
            shuffle!(aln, dim, true)
            @test aln != msa
            @test (aln .== GAP) == gaps[1]

            aln = copy(msa)
            shuffle!(aln, dim, false)
            @test aln != msa
            @test (aln .== GAP) != gaps[1]
        end

        @test_throws AssertionError shuffle(msa, 0, true)
        @test_throws AssertionError shuffle(msa, 3, true)
    end

    @testset "Columns" begin

        for i in 1:N
            # Fixed gaps
            msa = msas[i]
            aln = shuffle(msa, 2, true)
            @test aln != getresidues(msa)
            @test (aln .== GAP) == gaps[i]
            @test lcol[i] == mean(aln .== Residue('L'), dims=1)
            @test lseq[i] != mean(aln .== Residue('L'), dims=2)
            # Change gap positions
            aln = shuffle(msa, 2, false)
            @test aln != getresidues(msa)
            @test (aln .== GAP) != gaps[i]
            @test lcol[i] == mean(aln .== Residue('L'), dims=1)
        end
    end

    @testset "Sequences" begin

        for  i in 1:N
            # Fixed gaps
            msa = msas[i]
            aln = shuffle(msa, 1, true)
            @test aln != getresidues(msa)
            @test (aln .== GAP) == gaps[i]
            @test lcol[i] != mean(aln .== Residue('L'), dims=1)
            @test lseq[i] == mean(aln .== Residue('L'), dims=2)
            # Change gap positions
            aln = shuffle(msa, 1, false)
            @test aln != getresidues(msa)
            @test (aln .== GAP) != gaps[i]
            @test lseq[i] == mean(aln .== Residue('L'), dims=2)
        end
    end
end
