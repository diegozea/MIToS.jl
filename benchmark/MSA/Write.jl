let
    msa_file = joinpath(@__DIR__, "..", "..", "test", "data", "Gaoetal2011.fasta")
    msa = read_file(msa_file, FASTA, MultipleSequenceAlignment)

    SUITE["MSA"]["write"]["FASTA"] = @benchmarkable begin
        (tmp, io) = mktemp()
        close(io)
        write_file(tmp, $msa, FASTA)
        rm(tmp, force = true)
    end
end
