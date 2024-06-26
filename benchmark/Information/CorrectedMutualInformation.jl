let msa = rand(Random.MersenneTwister(1), res"ARNDCQEGHILKMFPSTWYV-", 50, 50),
    msa_large = msa[:, 1:10],
    msa_wide = msa[1:10, :]

    SUITE["Information"]["CorrectedMutualInformation"]["buslje09"]["msa"] = @benchmarkable buslje09($msa)
    SUITE["Information"]["CorrectedMutualInformation"]["buslje09"]["msa_large"] = @benchmarkable buslje09($msa_large)
    SUITE["Information"]["CorrectedMutualInformation"]["buslje09"]["msa_wide"] = @benchmarkable buslje09($msa_wide)
end