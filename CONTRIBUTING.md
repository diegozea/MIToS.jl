# Contributing

MIToS is a **Open Source** project, and you can contribute to it in different ways.
Please use [GitHub issues](https://github.com/diegozea/MIToS.jl/issues) to
**report errors/bugs** or to **ask for new features**. We welcome contributions in the
form of **pull requests**. For your code to be considered it must meet the
following guidelines.

  - By making a pull request, you agree to license your code under an MIT license.
  - Types and functions must be documented using Julia's docstrings.
  - All significant codes must be tested.

## Style

  - Type names are camel case, with the first letter capitalized.
    E.g. `MultipleSequenceAlignment`.
  - Function names, apart from constructors, are all lowercase. Include underscores between
    words only if the name would be hard to read without. E.g. `frequencies`, `read_file`.
  - Names of private (unexported) functions begin with an underscore, for example
    `_load_sequences`.
  - Separate logical blocks of code with blank lines.

### Code

MIToS has a `.JuliaFormatter.toml` file, so that [JuliaFormatter]() can be used to
automatically format the code following the described style.

  - Generally, keep lines below 92 columns.
  - Try to use a 4 spaces indentation.

### Documentation

  - Please include examples or `jldoctest` blocks if possible.

### References

Please include references to the papers where the algorithms are described. MIToS uses
[DocumenterCitations](https://github.com/JuliaDocs/DocumenterCitations.jl) to include
references in the documentation. All the references are stored in the `docs/src/refs.bib`
using the *BibTeX* format. When storing a new reference to the `refs.bib` file:

  - Please include the DOI of the paper so that the reference can have a link to the paper.
  - Use the paper's DOI as the reference key (*citekey*).

**To include a reference in the documentation**, use the `@cite` or the `@citet` macro.
The first will be rendered as a number, and the second as the author's name and the number.
For example, to include the reference to the MIToS paper as `[1]`, use
`[10.1093/bioinformatics/btw646](@cite)` in the documentation—please note that the DOI is
the citekey. To include it as `Zea et al. [1]`, use
`[10.1093/bioinformatics/btw646](@citet)`.

**If the reference is placed in a docstring**, to ensure that the reference is well rendered
in the REPL, please add it by hand using the first author's last name. Format it
using italic, for example, `*Zea et al.*`. Then add at the end of the docstring the
`# References` header, followed by a list of the references in MLA format. The whole
reference should link to the references section in the documentation.
Use `[MLA](@cite DOI)` to achieve that. For example:

```
# References

  - [Zea, Diego J., et al. "MIToS. jl: mutual information tools for protein sequence 
    analysis in the Julia language." Bioinformatics 33.4 (2017): 
    564-565.](@cite 10.1093/bioinformatics/btw646)
```

## Conduct

We adhere to the [Julia community standards](http://julialang.org/community/standards/).
