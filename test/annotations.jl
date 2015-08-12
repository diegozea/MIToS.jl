using Base.Test
using MIToS.MSA

print("""

Tests for Annotations
=====================
""")

print("""

Annotations() for empty annotations
""")
const annot = Annotations()
@test length(annot.file) == 0
@test length(annot.sequences) == 0
@test length(annot.columns) == 0
@test length(annot.residues) == 0

print("""
ncolums
""")
@test ncolumns(annot) == -1

print("""

Getters & Setters
-----------------
""")
setannotfile!(annot, "AC", "PF00571")
setannotcolumn!(annot, "SS_cons", "CCCCCHHHHHHHHHHHHHEEEEEEEEEEEEEEEEEEH")
setannotsequence!(annot, "O31698/88-139", "OS", "Bacillus subtilis")
setannotresidue!(annot, "O31698/18-71", "SS", "CCCHHHHHHHHHHHHHHHEEEEEEEEEEEEEEEEHHH")

@test getannotfile(annot, "AC") == "PF00571"
@test getannotcolumn(annot, "SS_cons") == "CCCCCHHHHHHHHHHHHHEEEEEEEEEEEEEEEEEEH"
@test getannotsequence(annot, "O31698/88-139", "OS") == "Bacillus subtilis"
@test getannotresidue(annot, "O31698/18-71", "SS") == "CCCHHHHHHHHHHHHHHHEEEEEEEEEEEEEEEEHHH"

@test ncolumns(annot) == 37
@test_throws DimensionMismatch setannotresidue!(annot, "O31698/18-71", "AS", "__*__")

print("""

copy, empty, isempty & empty!
-----------------------------
""")
copy_annot = copy(annot)
empty!(copy_annot)
@test ncolumns(annot) == 37
@test ncolumns(copy_annot) == -1

@test !isempty(annot)
@test isempty(copy_annot)

no_annotations = empty(Annotations)
@test ncolumns(no_annotations) == -1
@test isempty(no_annotations)

print("""

Filters
-------
""")
filtersequences!(annot, IndexedVector(["O31698/88-139", "O31698/18-71"]), [false, true])
@test length( getannotsequence(annot) ) == 0

mask = collect("CCCCCHHHHHHHHHHHHHEEEEEEEEEEEEEEEEEEH") .!= 'E'
filtercolumns!(annot, mask)
@test ncolumns(annot) == 19
@test getannotcolumn(annot, "SS_cons") == "CCCCCHHHHHHHHHHHHHH"

filtercolumns!(annot, [1,2,19])
@test ncolumns(annot) == 3
@test getannotcolumn(annot, "SS_cons") == "CCH"
