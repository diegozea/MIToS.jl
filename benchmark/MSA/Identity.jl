let
    msa_file = joinpath(@__DIR__, "..", "..", "test", "data", "Gaoetal2011.fasta")
    msa = read_file(msa_file, FASTA, MultipleSequenceAlignment)

    SUITE["MSA"]["identity"]["matrix_Float64"] =
        @benchmarkable percentidentity($msa, Float64)
    SUITE["MSA"]["identity"]["mean"] = @benchmarkable meanpercentidentity($msa)
end
