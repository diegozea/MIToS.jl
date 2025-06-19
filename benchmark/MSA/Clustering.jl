let
    msa_file = joinpath(@__DIR__, "..", "..", "test", "data", "Gaoetal2011.fasta")
    msa = read_file(msa_file, FASTA, MultipleSequenceAlignment)

    SUITE["MSA"]["hobohmI"]["pid62"] = @benchmarkable hobohmI($msa, 62)
end
