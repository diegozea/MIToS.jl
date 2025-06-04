# Contributor Guide

You will find more information about code style in the following file:
- /home/diego.zea/.julia/dev/MIToS/CONTRIBUTING.md

# Testing Instructions
To run the tests for this package, you can use the following command:

```bash
# From the repository root
julia --project -e 'using Pkg; Pkg.test(coverage=true)'
```

However, that run all the tests in the repository, which can take a long time. If you want 
to run only the tests for a specific module, `MSA` for example, you can use the following command:

```bash
julia --project -e 'push!(LOAD_PATH, "test"); using MIToSTests; MIToSTests.retest("MSA"); MIToSTests.retest("MSA")'
```
Note that the `MIToSTests.retest` function should be run two times. You can inspect the
output to look for the number of successful tests and the number of failed tests.