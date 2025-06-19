let
    msa_file = joinpath(@__DIR__, "..", "..", "test", "data", "Gaoetal2011.fasta")
    msa = read_file(msa_file, FASTA)

    SUITE["Information"]["highlevel"]["buslje09"] = @benchmarkable buslje09($msa)
    SUITE["Information"]["highlevel"]["BLMI"] = @benchmarkable BLMI($msa)
end
