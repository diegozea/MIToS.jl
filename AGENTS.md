# Contributor Guide

You will find more information about code style in the `CONTRIBUTING.md` file.
Follow the existing naming and file organization. Keep private helper names prefixed
with `_`, and group exports under the comment naming the file that defines them.

# Implementation

- Keep changes focused and small. Reuse existing helpers and public interfaces, and
  introduce abstractions only when they simplify the implementation.
- Prefer multiple dispatch over runtime inspection of method tables. Support built-in
  and user-defined types through the same documented interfaces.
- Use authoritative specifications and official documentation to establish expected
  behavior. Write original implementations; do not copy code from external sources,
  regardless of project or license.

# Workflow

- Preserve existing work, including uncommitted changes. Do not discard, overwrite or
  commit the user's changes without authorization.
- Respect the requested boundaries between planning, local changes, publishing and
  merging. Iterate locally when requested, and merge only with explicit authorization.
- When verifying CI, check the latest relevant commit and all applicable checks,
  including individual coverage reports. After merging, verify the checks on the
  resulting commit on the target branch.

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

Run `MIToSTests.retest` twice and verify that the selected tests actually ran: zero tests
is not a successful validation. Report the numbers of passing, failing and skipped tests.

- Follow the existing test structure. Keep helpers inside the relevant `@testset` blocks,
  and avoid new global test definitions or additional test modules.
- Centralize `using` and `import` statements in `test/tests.jl`. Add short comments naming
  the test files that need those imports, without repeating the directory path.
- Briefly explain the behavior each testset protects. For regressions, reproduce the
  failure before fixing it when practical, and compare against the PR's base branch.
- Keep tests focused on behavior. Merge redundant tests without losing coverage, and
  add compact tests for documentation examples when they protect meaningful behavior.
- Give synthetic fixtures names that cannot be mistaken for real datasets.

Run tests that download data with network access. If the network or a required service is
unavailable, skip only the affected tests and report the skips. Do not hide parsing or
other functional failures as connectivity problems. Avoid adding servers, socket handling
or other test infrastructure merely to make network tests run offline.

# Dependencies

Install new dependencies in the appropriate environment before running tests. Keep runtime,
test and development-tool dependencies separate, following the repository's existing
dependency configuration. Make tools such as ReTest, JuliaFormatter and benchmarking tools
available from a development environment; do not add them as runtime dependencies merely
to prepare local checks.

# Formatting

Run JuliaFormatter on changed Julia files, using the repository's configuration. Do not
approximate its formatting manually. For example, to format `abc.jl`:

```bash
julia --project -e 'using JuliaFormatter; JuliaFormatter.format_file("abc.jl")'
```

If formatting changes docstring escaping or layout, verify the rendered documentation.

# Documentation

- Write user guides for computational biologists, including readers coming from biology
  with limited programming experience. Keep explanations concise and avoid unnecessary
  computing jargon or references to a PR's history.
- Show complete recommended examples and omit redundant default keyword arguments.
  Keep detailed API contracts in docstrings and extension instructions in the Development
  documentation.
- Keep warnings and error messages short and informative.

# Release Notes

Edit `NEWS.md` or change the version in `Project.toml` only when explicitly requested.
The package follows semantic versioning. Each section in `NEWS.md` has a title indicating the
previous and current versions. The last version should always be placed first in the 
`NEWS.md` file, followed by older sections, ordered from most recent to oldest. Document 
each change with bullet points. Clearly label breaking changes using 
the `*[Breaking change]*` tag at the beginning of the bullet, so they are easily
identifiable.

# Benchmarking

Base performance claims on measurements against the relevant baseline, normally the PR's
base branch. Inspect timings and allocations rather than relying only on a green CI job.

If you are explicitly asked to run the benchmark suite, make sure
`PkgBenchmark` and `BenchmarkTools` are installed. Then execute the
following command from the repository root to tune and run all benchmarks:

```bash
julia --project -e 'import PkgBenchmark, MIToS; PkgBenchmark.benchmarkpkg(MIToS; retune=true)'
```

This command creates a `benchmark/tune.json` file with the tuning
information and prints benchmark results to the terminal.
