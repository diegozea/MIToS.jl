# Run from the repository root with MIToS and BenchmarkTools on LOAD_PATH:
# julia --project --threads=4,1 benchmark/threaded_mapcolpairfreq.jl
using BenchmarkTools
using MIToS.MSA
using MIToS.Information
using MIToS.Utils
using Random
using Statistics

function benchmark_pairs(label, msa)
    serial_table = Frequencies(ContingencyTable(Float64, Val{2}, UngappedAlphabet()))
    threaded_table = deepcopy(serial_table)
    serial() = mapcolpairfreq!(
        normalized_mutual_information,
        msa,
        serial_table;
        usediagonal = false,
    )
    threaded() = mapcolpairfreq!(
        normalized_mutual_information,
        msa,
        threaded_table;
        usediagonal = false,
        threads = true,
    )
    @assert isequal(serial(), threaded())
    # Include output/scratch allocation, but exclude compilation and input construction.
    serial_trial = @benchmark $serial() samples = 30 seconds = 3 evals = 1
    threaded_trial = @benchmark $threaded() samples = 30 seconds = 3 evals = 1
    s, t = median(serial_trial), median(threaded_trial)
    println(
        join(
            (
                label,
                size(msa, 1),
                size(msa, 2),
                s.time / 1e6,
                t.time / 1e6,
                s.time / t.time,
                s.memory,
                t.memory,
            ),
            '\t',
        ),
    )
end

println(
    "Julia ",
    VERSION,
    "; worker threads=",
    Threads.nthreads(:default),
    "; interactive threads=",
    Threads.nthreads(:interactive),
    "; CPU=",
    Sys.CPU_NAME,
)
println(
    "case\tsequences\tcolumns\tserial_ms\tthreaded_ms\tserial/threaded\tserial_bytes\tthreaded_bytes",
)
data = joinpath(@__DIR__, "..", "test", "data")
benchmark_pairs("Gaoetal2011", read_file(joinpath(data, "Gaoetal2011.fasta"), FASTA))
benchmark_pairs("PF09645", read_file(joinpath(data, "PF09645_full.fasta.gz"), FASTA))
for (label, nseq, ncol) in (
    ("synthetic_medium", 500, 200),
    ("synthetic_deep", 5000, 100),
    ("synthetic_wide", 100, 1000),
)
    msa = rand(MersenneTwister(197), res"ARNDCQEGHILKMFPSTWYV-", nseq, ncol)
    benchmark_pairs(label, msa)
end
