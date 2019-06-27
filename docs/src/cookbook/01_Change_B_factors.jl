# # Change B-factors
#
#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__cookbook/notebooks/01_Change_B_factors.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__cookbook/notebooks/01_Change_B_factors.ipynb)
#
#
# ## Problem description
#
# It is a common practice to change the B-factors of a PDB to store information
# about atoms or residues to be used by other programs. In particular, values
# in the B-factor column can be easily used to colour residues with
# [PyMOL](https://pymolwiki.org/index.php/Color#B-Factors) or
# [Chimera](https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/bfactor.html).
#
# We cannot simply assign a new value to the `B` field of a `PDBAtom` because
# this type is immutable. However, we can make use of the `@set` macro of the
# [Setfield](https://github.com/jw3126/Setfield.jl) package to create a new
# `PDBAtom` with a different B-factor value.
#
# In a PDB file, B-factors are stored from the column 61 to 66. Therefore, new
# B-factors should be a `String` with 6 or fewer characters, normally using two
# characters for decimal values. We can use `fmt` and `FormatSpec` from the
# [Formatting](https://github.com/JuliaIO/Formatting.jl) package to create a
# proper B-factor string.
#
# ## MIToS solution
#
# For this example we are going to use the small heat shock protein AgsA from
# *Salmonella typhimurium* (PDB code: *4ZJ9*) available in MIToS docs data:

using MIToS
pdbfile = abspath(pathof(MIToS), "..", "..", "docs", "data", "4zj9.pdb")

# First, we need to read the PDB file using the `MIToS.PDB` module:

using MIToS.PDB
residues = read(pdbfile, PDBFile)

# ...

# ## Further details
