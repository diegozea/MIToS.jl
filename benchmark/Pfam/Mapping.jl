let pfam_alignment = read_file(
        joinpath(@__DIR__, "..", "..", "docs", "data", "PF18883.stockholm.gz"),
        Stockholm,
    )
    SUITE["Pfam"]["accession mapping"]["acc2seqnames"] =
        @benchmarkable MIToS.Pfam._get_acc2seqnames($pfam_alignment)
end
