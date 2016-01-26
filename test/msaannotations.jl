# using Base.Test
# using MIToS.MSA

print("""
Tests for Mappings on Annotated Multiple Sequence Alignments
------------------------------------------------------------
""")

let pfam = read(joinpath(pwd(), "data", "PF09645_full.stockholm"), Stockholm,
                generatemapping=true, useidcoordinates=true),
    F112_SSV1 = collect(string(".....QTLNSYKMAEIMYKILEKKGELTLEDILAQFEISVPSAYNIQRALKAICERHPDECEVQY",
                               "KNRKTTFKWIKQEQKEEQKQEQTQDNIAKIFDAQPANFEQTDQGFIKAKQ....."))

  @test minimum(getsequencemapping(pfam, 4)) == 3
  @test maximum(getsequencemapping(pfam, "F112_SSV1/3-112")) == 112
  @test getcolumnmapping(pfam) == getindex([ i for i in 1:length(F112_SSV1) ], F112_SSV1 .!= '.')
  @test length(pfam.annotations.sequences) == 9
end

let fasta = read(joinpath(pwd(), "data", "Gaoetal2011.fasta"), FASTA, generatemapping=true)

  @test getsequencemapping(fasta, 1) == [1, 2, 3, 4, 5, 6]
  @test getsequencemapping(fasta, "SEQ2") == [1, 2, 3, 4, 5, 6]
  @test getcolumnmapping(fasta) == [1, 2, 3, 4, 5, 6]
end

print("""

Tests for annotation of modifications
-------------------------------------
""")

let pfam = read(joinpath(pwd(), "data", "PF09645_full.stockholm"), Stockholm,
                generatemapping=true, useidcoordinates=true), io = IOBuffer()

  print("""
  Test printing annotations and modifications
  """)
  print(io, pfam.annotations)
  printed = split(takebuf_string(io), '\n')
  @test length(printed) == 17
  @test printed[1] == string("#=GF NCol   120")
  @test printed[2] == string("#=GF ColMap   6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22",
                             ",23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,",
                             "43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,",
                             "63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,",
                             "84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,",
                             "103,104,105,106,107,108,109,110,111,112,113,114,115")
  @test ismatch(r"MIToS_", printed[3])
  @test ismatch(r"MIToS_", printed[4])

  print("""
  Test delete MIToS annotations
  """)
  delete_annotated_modifications!(pfam)
  print(io, pfam.annotations)
  printed = split(takebuf_string(io), '\n')
  @test length(printed) == 17 - 2 # 2 MIToS annotations
  @test !ismatch(r"MIToS_", printed[3])
  @test !ismatch(r"MIToS_", printed[4])
end
