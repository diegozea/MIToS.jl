# Contributor Guide

You will find more information about code style in the following file:

* /home/diego.zea/.julia/dev/MIToS/CONTRIBUTING.md

# Testing Instructions

To run the tests for this package, you can use the following command:

```bash
# From the repository root
julia --project -e 'using Pkg; Pkg.test(coverage=true)'
```

However, that runs all the tests in the repository, which can take a long time. If you want
to run a specific `@testset` named `abc`, for example, you can use the following command:

```bash
julia --project -e 'push!(LOAD_PATH, "test"); using MIToSTests; MIToSTests.retest("abc"); MIToSTests.retest("abc")'
```

Note that the `MIToSTests.retest` function should be run two times. You can inspect the
output to look for the number of successful tests and the number of failed tests. That is
important as some tests do not fail because they don't run (the count of tests is 0).

# Formatting

At the end, you can format your files, e.g., the `abc.jl` file, using `JuliaFormatter`:

```bash
julia --project -e 'using JuliaFormatter; JuliaFormatter.format_file("abc.jl")'
```
