let
    msa_file = joinpath(@__DIR__, "..", "..", "test", "data", "Gaoetal2011.fasta")
    msa = read_file(msa_file, FASTA, MultipleSequenceAlignment)

    SUITE["MSA"]["hobohmI"]["pid20"] = @benchmarkable hobohmI($msa, 20)
    SUITE["MSA"]["hobohmI"]["pid62"] = @benchmarkable hobohmI($msa, 62)
    SUITE["MSA"]["hobohmI"]["pid80"] = @benchmarkable hobohmI($msa, 80)
    SUITE["MSA"]["hobohmI"]["pid99"] = @benchmarkable hobohmI($msa, 99)
end
