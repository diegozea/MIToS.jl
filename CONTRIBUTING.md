# Contributing

MIToS is and **Open Source** project and there are different ways to contribute.
Please, use [GitHub issues](https://github.com/diegozea/MIToS.jl/issues) to
**report errors/bugs** or to **ask for new features**.
We welcome contributions in the form of **pull requests**. For your code to be considered
it must meet the following guidelines.  

- By making a pull request, you're agreeing to license your code under a MIT license.
- Types and functions must be documented using Julia's docstrings.
- All significant code must be tested.

## Style

- Type names are camel case, with the first letter capitalized.
E.g. `MultipleSequenceAlignment`.
- Function names, apart from constructors, are all lowercase. Include underscores between
words only if the name would be hard to read without. E.g. `covalentradius`,
 `check_atoms_for_interactions`.
- Names of private (unexported) functions begin with underscore.
- Separate logical blocks of code with blank lines.
- Generally try to keep lines below 92-columns, unless splitting a long line onto multiple
lines makes it harder to read.
- Try to use a 4 spaces indentation

### Documentation

- Do not include headers/titles in the docstrings
- Please include examples if it's possible

## Conduct

We adhere to the [Julia community standards](http://julialang.org/community/standards/).
