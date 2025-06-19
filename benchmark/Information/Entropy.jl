let
    data_dir = joinpath(@__DIR__, "..", "..", "test", "data")
    msa_file = joinpath(data_dir, "PF09645_full.fasta.gz")
    msa = read_file(msa_file, FASTA)
    table = Probabilities{Float64,1,UngappedAlphabet}(
        ContingencyTable(Float64, Val{1}, UngappedAlphabet()),
    )

    SUITE["Information"]["shannon_entropy"]["PF09645"] =
        @benchmarkable mapcolfreq!(shannon_entropy, $msa, $table)
end
