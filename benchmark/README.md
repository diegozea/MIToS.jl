To run the benchmark suite, you must have `PkgBenchmark` and `BenchmarkTools` installed. 
Then, you can run the following code in the Julia REPL:

```julia
import PkgBenchmark, BenchmarkTools, MIToS; PkgBenchmark.benchmarkpkg(MIToS)

```

For a direct serial/threaded comparison of `mapcolpairfreq!`, run
`julia --project --threads=4,1 benchmark/threaded_mapcolpairfreq.jl` with
`BenchmarkTools` available in your environment. The script checks exact score equality
and reports median times and allocations for two real fixtures and synthetic alignments
with different sequence and column counts. Input generation and compilation are excluded;
output and scratch-table allocations are included. Run without competing CPU workloads.
Repeat with `--threads=1` to check the serial fallback. Results depend on alignment size,
callback cost, hardware and available CPU time; small inputs may be slower with threading.
