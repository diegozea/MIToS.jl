# Contributor Guide

You will find more information about code style in the `CONTRIBUTING.md` file.

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

If your modifications introduce a new dependency, you should install it before running
the tests. You can do that by running the following command:

```bash
julia --project -e 'using Pkg; Pkg.add("NewDependency")'
```

# Formatting

At the end, you can format your files, e.g., the `abc.jl` file, using `JuliaFormatter`:

```bash
julia --project -e 'using JuliaFormatter; JuliaFormatter.format_file("abc.jl")'
```

# Release Notes

The `NEWS.md` file contains the release notes for this package. Please update it with any changes you make to the package. Note that this package follows semantic versioning, so any changes that break the API should increment the major version number with respect to the current version. You can check the current version in the `Project.toml` file. Please, do not update the version number in the `Project.toml` file.

Each new set of changes should be documented under a section titled:

```md
### Changes from CURRENT to master
```

For example, if the current version listed in the `Project.toml` file is `0.1.0`, the title of your next set of changes should be:

```md
### Changes from 0.1.0 to master
```

This section should always be placed first in the `NEWS.md` file, followed by older sections, ordered from most recent to oldest.

Document each change with bullet points. Clearly label breaking changes using the `*[Breaking change]*` tag at the beginning of the bullet, so they are easily identifiable.
