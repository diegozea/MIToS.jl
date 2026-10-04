# AlphaFold mmCIF and template-hit-name checks

`AlphaFoldMMCIF.jl` exports four templates and checks them with AlphaFold's actual
`mmcif_parsing.parse`, including the sequence-to-coordinate mapping. It also passes
ColabFold-style hit names to AlphaFold's actual `_get_pdb_id_and_chain` and verifies
that the returned identifier locates the exported file. It covers a PDB
with insertion codes, explicit author renumbering, model 14 of an NMR ensemble, and
multiple chains with a modified HETATM residue. It runs separately from `Pkg.test()`
and does not add runtime or test dependencies to MIToS. The filenames (`inst.cif`,
`tmpl.cif`, `nmr1.cif`, `mods.cif`) use four-character lowercase alphanumeric stems.
Their descriptive mmCIF `entry_id` values are independent of those filenames.

With an AlphaFold checkout and its parser dependencies available to Python:

```sh
PYTHONPATH=/path/to/alphafold \
ALPHAFOLD_PYTHON=/path/to/python \
julia --project test/integration/AlphaFoldMMCIF.jl
```

The check was exercised on Linux with Julia 1.13.1, Biopython 1.88, and AlphaFold
commit `c77e5d2a8961d1a353632c462914ff0a32a950f6`. That parser imports `absl-py`,
`biopython`, `numpy` and `jax`; model weights are unnecessary. Warnings about missing
release dates are expected: the export must not invent those dates.

AlphaFold inference, HHsearch, Kalign and complete ColabFold execution remain outside
this check. ColabFold's insertion handling, fallback dates and
modified-residue handling are discussed in the PDB documentation.
