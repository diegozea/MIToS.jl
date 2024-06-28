let seq_a = rand(Random.MersenneTwister(37), res"ARNDCQEGHILKMFPSTWYV-", 500),
    seq_b = rand(Random.MersenneTwister(73), res"ARNDCQEGHILKMFPSTWYV-", 500),
    Na = ContingencyTable(Float64, Val{1}, UngappedAlphabet()),
    Nab = ContingencyTable(Float64, Val{2}, UngappedAlphabet())

    SUITE["Information"]["frequencies!"]["1"] =
        @benchmarkable frequencies!($Na, $seq_a)
    SUITE["Information"]["frequencies!"]["2"] =
        @benchmarkable frequencies!($Nab, $seq_a, $seq_b)
end
