let
    sifts_file = joinpath(@__DIR__, "..", "..", "test", "data", "2vqc.xml.gz")

    SUITE["SIFTS"]["siftsmapping"]["2vqc"] = @benchmarkable SIFTS.siftsmapping(
        $sifts_file,
        SIFTS.dbPDBe,
        "2vqc",
        SIFTS.dbPDB,
        "2vqc";
        chain = "A",
        missings = false,
    )
end
