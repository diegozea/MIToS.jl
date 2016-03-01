# Contributing
MIToS is and **Open Source** project and there are different ways to contribute.
Please, use [GitHub issues](https://github.com/diegozea/MIToS.jl/issues) to **report bugs** or to **ask for new features**.
We welcome contributions in the form of **pull requests**. For your code to be considered it must meet the following guidelines.
- By making a pull request, you're agreeing to license your code under a MIT license.
- Types and functions must be documented using Julia's docstrings.
- All significant code must be tested.

## Style
- Indent with 2 spaces.
- Type names are camel case, with the first letter capitalized. E.g. `MultipleSequenceAlignment`.
- Function names, apart from constructors, are all lowercase. Include underscores between words only if the name would be hard to read without. E.g. `covalentradius`, `check_atoms_for_interactions`.
- Generally try to keep lines below 80-columns, unless splitting a long line onto multiple lines makes it harder to read.
- Separate logical blocks of code with blank lines.
- Document functions using bare docstrings before a definition:
- Functions that get or set variables in a type should not be prefixed with 'get' or 'set'. The getter should be named for the variable it sets, and the setter should have the same name as the getter, with the suffix !.

## Conduct
We adhere to the [Julia community standards](http://julialang.org/community/standards/).
