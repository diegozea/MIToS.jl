# using Base.Test
# using MIToS.MSA

print("""

Tests for Annotations
=====================
""")

print("""

Annotations() for empty annotations
""")
let annot = Annotations()
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
  setannotresidue!(annot, "O31698/18-71", "SS", "CCCHHHHHHHHHHHHHHHEEEEEEEEEEEEEEEEHHH")

  @test_throws ErrorException setannotresidue!(annot, "O31698/18-71", randstring(51),    "CCCHHHHHHHHHHHHHHHEEEEEEEEEEEEEEEEHHH")
  @test_throws ErrorException setannotresidue!(annot, "O31698/18-71", "My Feature Name", "CCCHHHHHHHHHHHHHHHEEEEEEEEEEEEEEEEHHH")

  @test ncolumns(annot) == 37

  setannotfile!(annot, "AC", "PF00571")
  setannotcolumn!(annot, "SS_cons", "CCCCCHHHHHHHHHHHHHEEEEEEEEEEEEEEEEEEH")
  setannotsequence!(annot, "O31698/88-139", "OS", "Bacillus subtilis")

  @test getannotfile(annot, "AC") == "PF00571"
  @test getannotcolumn(annot, "SS_cons") == "CCCCCHHHHHHHHHHHHHEEEEEEEEEEEEEEEEEEH"
  @test getannotsequence(annot, "O31698/88-139", "OS") == "Bacillus subtilis"
  @test getannotresidue(annot, "O31698/18-71", "SS") == "CCCHHHHHHHHHHHHHHHEEEEEEEEEEEEEEEEHHH"

  @test getannotfile(annot, "An", "Default") == "Default"
  @test getannotcolumn(annot, "Other", "Default") == "Default"
  @test getannotsequence(annot, "O31698/1-88", "OS", "Default") == "Default"
  @test getannotresidue(annot, "O31698/1-88", "SS", "Default") == "Default"

  @test ncolumns(annot) == 37
  @test_throws DimensionMismatch setannotresidue!(annot, "O31698/18-71", "AS", "__*__")
  @test_throws DimensionMismatch setannotcolumn!(annot, "SS_cons", "---CCCCCHHHHHHHHHHHHHEEEEEEEEEEEEEEEEEEH---")

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
  filtersequences!(annot, IndexedArray(["O31698/88-139", "O31698/18-71"]), [false, true])
  @test length( getannotsequence(annot) ) == 0
  filtersequences!(annot, IndexedArray(["O31698/88-139", "O31698/18-71"]), [true, false])
  @test length( getannotresidue(annot) ) == 0

  mask = collect("CCCCCHHHHHHHHHHHHHEEEEEEEEEEEEEEEEEEH") .!= 'E'
  filtercolumns!(annot, mask)
  @test ncolumns(annot) == 19
  @test getannotcolumn(annot, "SS_cons") == "CCCCCHHHHHHHHHHHHHH"

  filtercolumns!(annot, [1,2,19])
  @test ncolumns(annot) == 3
  @test getannotcolumn(annot, "SS_cons") == "CCH"

end
