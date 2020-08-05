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
#md nothing # hide

# First, we need to read the PDB file using the `MIToS.PDB` module:

using MIToS.PDB
pdb_residues = read(pdbfile, PDBFile)
#md nothing # hide

# For this example, we are going to replace the B-factor of the alpha-carbons
# by the residue hydrophobicity according to the hydrophobicity scale of
# [Kyte and Doolittle](https://doi.org/10.1016/0022-2836(82)90515-0) used by
# [Chimera](https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/midas/hydrophob.html):

hydrophobicity = Dict(
	"ILE" => 4.5,
	"VAL" => 4.2,
	"LEU" => 3.8,
	"PHE" => 2.8,
	"CYS" => 2.5,
	"MET" => 1.9,
	"ALA" => 1.8,
	"GLY" => -0.4,
	"THR" => -0.7,
	"SER" => -0.8,
	"TRP" => -0.9,
	"TYR" => -1.3,
	"PRO" => -1.6,
	"HIS" => -3.2,
	"GLU" => -3.5,
	"GLN" => -3.5,
	"ASP" => -3.5,
	"ASN" => -3.5,
	"LYS" => -3.9,
	"ARG" => -4.5 )
#md nothing # hide

# First, we define a helper function using `Formatting` to create a proper
# B-factor string with the PDB format; 6 characters and 2 digits after the
# decimal point.
# The [PDB format description](https://www.wwpdb.org/documentation/file-format-content/format23/sect9.html)
# describe this field as:
# ```
# COLUMNS      DATA TYPE        FIELD      DEFINITION
# ------------------------------------------------------
# 61 - 66      Real(6.2)        tempFactor Temperature factor.
# ```

using Formatting

"Return value as a string with the B factor format described in PDB."
format_b_factor(value) = fmt(FormatSpec("6.2f"), value) # e.g. 1.5 -> "  1.50"
#md nothing # hide

# Then, where are using that helper function to define a function that returns
# a new `PDBAtom` by changing the `B` factor field using the `Setfield` package.

using Setfield

"Return a new PDBAtom with the B-factor changed to value."
function change_b_factor(atom::PDBAtom, value)
	b_factor_string = format_b_factor(value)
	b_factor_string = strip(b_factor_string) # e.g. "  1.50" -> "1.50"
	if length(b_factor_string) > 6
		throw(ErrorException("$b_factor_string has more than 6 characters."))
	end
	@set atom.B = b_factor_string
end
#md nothing # hide

# Now, we can use the `change_b_factor` function to change the B-factor of each
# `"CA"` atom:

for res in pdb_residues
	for i in eachindex(res.atoms)
		atom = res.atoms[i]
		if atom.atom == "CA"
			res.atoms[i] = change_b_factor(atom, hydrophobicity[res.id.name])
		end
	end
end

# Finally, we can save the changed residues in a new PDB file.
#
# ```julia
# write("4zj9_hydrophobicity.pdb", pdb_residues, PDBFile)
# ```
#
# ## Discussion
#
# While we have focused on changing the B-factor field of a `PDBAtom`, you can
# use the same approach to change other fields. However, if you want to change
# atom coordinates, it is better to use the `change_coordinates` function from
# the PDB module of MIToS.
#
# MIToS atoms and residues generally stores the string present in the input
# file without surrounding spaces. You can use the `Formatting` module to
# create these strings and `strip` to get rid of the spaces. You can see the
# [PDB format description](https://www.wwpdb.org/documentation/file-format-content/format23/sect9.html)
# to know what is the format of the expected string or see the
# [MIToS PDB print source code](https://github.com/diegozea/MIToS.jl/blob/master/src/PDB/PDBParser.j)
# to get a quick idea.
