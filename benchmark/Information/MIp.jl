let
    data_dir = joinpath(@__DIR__, "..", "..", "test", "data")
    msa_file = joinpath(data_dir, "PF09645_full.fasta.gz")
    msa = read_file(msa_file, FASTA, MultipleSequenceAlignment)
    table = Counts{Float64,2,GappedAlphabet}(
        ContingencyTable(Float64, Val{2}, GappedAlphabet()),
    )

    SUITE["Information"]["MIp"]["PF09645"] = @benchmarkable begin
        mi = mapcolpairfreq!(mutual_information, $msa, $table)
        APC!(mi)
    end
end
