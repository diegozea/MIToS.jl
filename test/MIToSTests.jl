"""
You need to `include` this file or load this module in your Julia session to run the tests
using `ReTest`. The main advantage of `ReTest` is that it allows you to run only selected
tests. As `ReTest` is not a `MIToS` dependency, so you need to install it manually.
It is also recommended to install `Revise`. Also, you will need to install the test
dependencies, `Documenter` and `ROCAnalysis`, outside the MIToS environment:

```julia
using Pkg
Pkg.add("ReTest")
Pkg.add("Revise")
Pkg.add("Documenter")
Pkg.add("ROCAnalysis")
```

An example of usage if you want to run the `hcat` tests from the `MSA` module:

```julia
push!(LOAD_PATH, joinpath(homedir(), ".julia", "dev", "MIToS", "test"))
using Revise, MIToSTests
MIToSTests.retest("MSA")
MIToSTests.retest("hcat")
```

Note that we need to fisrt run the most general test and then the specific one. Otherwise,
`ReTest` will not be able to find the `hcat` test.

NOTE: For some reason, after modifying the tests, `Revise` does not detect the changes
automatically. However, runing `retest` again for the whole module looks to do the trick.
"""
module MIToSTests
using ReTest
include("tests.jl")
end
