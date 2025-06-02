# Contributor Guide

You will find more information about code style in the following file:
- /home/diego.zea/.julia/dev/MIToS/CONTRIBUTING.md

# Testing Instructions
To run the tests for this package, you can use the following command:

```bash
# From the repository root
julia --project -e 'using Pkg; Pkg.test(coverage=true)'
```