using Base.Test
using MIToS.MSA

print("""

Tests for Mappings on Annotated Multiple Sequence Alignments
============================================================
""")

print("""

Read and generate mappings
--------------------------
""")

let pfam = read("./data/PF09645_full.stockholm", Stockholm,
                generatemapping=true, useidcoordinates=true),
    F112_SSV1 = collect(".....QTLNSYKMAEIMYKILEKKGELTLEDILAQFEISVPSAYNIQRALKAICERHPDECEVQYKNRKTTFKWIKQEQKEEQKQEQTQDNIAKIFDAQPANFEQTDQGFIKAKQ.....")

  @test minimum(getsequencemapping(pfam, 4)) == 3
  @test maximum(getsequencemapping(pfam, "F112_SSV1/3-112")) == 112
  @test getcolumnmapping(pfam) == getindex([ i for i in 1:length(F112_SSV1) ], F112_SSV1 .!= '.')
  @test length(pfam.annotations.sequences) == 9
end

let fasta = read("./data/Gaoetal2011.fasta", FASTA, generatemapping=true)

  @test getsequencemapping(fasta, 1) == [1, 2, 3, 4, 5, 6]
  @test getsequencemapping(fasta, "SEQ2") == [1, 2, 3, 4, 5, 6]
  @test getcolumnmapping(fasta) == [1, 2, 3, 4, 5, 6]
end
