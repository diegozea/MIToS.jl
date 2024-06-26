To run the benchmark suite, you need to have the `PkgBenchmark` package installed. 
Then, you can run the following code in the Julia REPL:

```julia
import PkgBenchmark
import MIToS
PkgBenchmark.benchmarkpkg(MIToS)
```