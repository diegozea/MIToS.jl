using Base.Test
using MIToS.MSA

# Fields of MultipleSequenceAlignment
const msa_fields = Symbol[:id, :msa, :sequencemapping, :filecolumnmapping, :annotations]
# Fields of AlignedSequence
const seq_fields = Symbol[:id, :index, :sequence, :sequencemapping, :filecolumnmapping, :annotations]

print("""

Tests for Multiple Sequence Alignments
======================================
""")

print("""

Parse Pfam
----------
""")

print("""
Test pfam stockholm parser using the 4 sequence full MSA for PF09645
> Order: Tree
> Inserts lower case
> Gaps as "." or "-" (mixed)
> Pfam version 28.0, based on UniProt release 2014_07
""")

const pfam = readpfam("./data/PF09645_full.stockholm")
const pfam_na = readpfam("./data/PF09645_full.stockholm", MultipleSequenceAlignment)
const F112_SSV1 = collect(".....QTLNSYKMAEIMYKILEKKGELTLEDILAQFEISVPSAYNIQRALKAICERHPDECEVQYKNRKTTFKWIKQEQKEEQKQEQTQDNIAKIFDAQPANFEQTDQGFIKAKQ.....")

@test pfam_na == convert(MultipleSequenceAlignment, pfam)

@test size(pfam.msa, 1) == 4
@test size(pfam.msa, 2)  == length(F112_SSV1[ F112_SSV1 .!= '.' ]) # Without inserts
@test slice(pfam.msa,4,:) == convert(Vector{Residue}, F112_SSV1[ F112_SSV1 .!= '.' ])
@test minimum(pfam.sequencemapping[4,:]) == 3
@test maximum(pfam.sequencemapping[4,:]) == 112
@test pfam.id.values == ["C3N734_SULIY/1-95", "H2C869_9CREN/7-104", "Y070_ATV/2-70", "F112_SSV1/3-112"]
@test pfam.filecolumnmapping.values == getindex([ i for i in 1:length(F112_SSV1) ], F112_SSV1 .!= '.')
@test !isempty(pfam.annotations)
@test length(pfam.annotations.file) == 0
@test length(pfam.annotations.sequences) == 5
@test length(pfam.annotations.columns) == 2
@test length(pfam.annotations.residues) == 1
@test pfam.annotations.residues[("F112_SSV1/3-112","SS")] == "X---HHHHHHHHHHHHHHHSEE-HHHHHHHH---HHHHHHHHHHHHHHHHH-TTTEEEEE-SS-EEEEE--XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
@test pfam.annotations.columns["seq_cons"] == "...NshphAclhaKILppKtElolEDIlAQFEISsosAYsI.+sL+hICEpH.-ECpsppKsRKTlhh.hKpEphppptpEp..ppItKIhsAp................"
@test pfam.annotations.sequences[("F112_SSV1/3-112","DR")] == "PDB; 2VQC A; 4-73;"

print("""

Parse Fasta
-----------
""")

const fasta = readfasta("./data/PF09645_full.fasta")
const small = readfasta("./data/Gaoetal2011.fasta")
const fasta_na = readfasta("./data/PF09645_full.fasta", MultipleSequenceAlignment)
const small_na = readfasta("./data/Gaoetal2011.fasta", MultipleSequenceAlignment, deletefullgaps=false)

@test fasta_na == convert(MultipleSequenceAlignment, fasta)
@test small_na == convert(MultipleSequenceAlignment, small)

print("""
Test the FASTA parser usando Gao et.al. 2011 example
""")

const seq1 = res"DAWAEE"
const seq2 = res"DAWAEF"
const seq3 = res"DAWAED"
const seq4 = res"DAYCMD"
const seq5 = res"DAYCMT"
const seq6 = res"DAYCMT"
const gaores = [seq1'; seq2'; seq3'; seq4'; seq5'; seq6']

@test small.msa == gaores
@test small.id.values == ASCIIString[ "SEQ$i" for i in 1:6 ]
@test small.sequencemapping == [1  2  3  4  5  6
                                1  2  3  4  5  6
                                1  2  3  4  5  6
                                1  2  3  4  5  6
                                1  2  3  4  5  6
                                1  2  3  4  5  6]
@test small.filecolumnmapping.values == [1, 2, 3, 4, 5, 6]
@test isempty(small.annotations)

print("""
Are the parsers (pfam/fasta) getting the same result?
""")

for field in msa_fields
  if field != :annotations
    @eval @test fasta.$field == pfam.$field
  end
end
@test isempty(fasta.annotations)

print("""

Selection without Mappings
--------------------------
""")

@test getresidues(small) == gaores
@test getresidues(small_na) == gaores

print("""
Test getindex
""")

for index in [ 2, (2,:), (:,2), (2,2), (2:4, 2:4), (2,[3,4,5]),
              (:, [1,2,3,4,5,6] .< 4 ), ([1,2,3,4,5,6] .< 4, :),
              (:,:), ([1,2,3,4,5,6] .< 4, [1,2,3,4,5,6] .< 4) ]
  @test small[index...] == gaores[index...]
  @test small_na[index...] == gaores[index...]
end

print("""
Test setindex! with copy and deepcopy
""")
copy_small = copy(small)
deepcopy_small_na = deepcopy(small_na)

for (index, value) in [ ((1),Residue('H')), # ((1),Residue('D')),
                        ((:,1), res"HHHHHH"), # ((:,1), res"DDDDDD"),
                        ((1),Residue('H')) ] #, ((1),Residue('D')) ]
  copy_small[index...] = value
  @test copy_small[index...] == value
  @test small[index...] != value
  deepcopy_small_na[index...] = value
  @test deepcopy_small_na[index...] == value
  @test small_na[index...] != value
end

for (index, value) in [ (1,Residue('H')), (1,Residue('D')),
                        ( : , res"HHHHHH"), ( : , res"DAYCMD") ]
  tmp_sequence = copy(getsequence(small, 4))
  tmp_sequence[index] = value
  @test tmp_sequence[index] == value
end

print("""
Test getresidues, getsequence and getresiduesequence
""")

for i in 1:6
  println("Test sequence $i for getresidues and getsequence")
  @eval @test getresidues(getsequence(small,$i)) == $(Symbol("seq$i"))
  @eval @test getresidues(getsequence(small_na,$i)) == $(Symbol("seq$i"))
end

const rawseqs = getresiduesequences(small)
const rawseqs_na = getresiduesequences(small_na)

for i in 1:6
  println("Test sequence $i from getresiduesequences")
  @eval @test rawseqs[$i] == $(Symbol("seq$i"))
  @eval @test rawseqs_na[$i] == $(Symbol("seq$i"))
end

print("""

Size, sequences and columns
---------------------------
""")

for aln in (small, small_na)
  @test size(aln) == (6,6)
  @test length(aln) == 36
  @test ncolumns(aln) == 6
  @test nsequences(aln) == 6
end

@test ncolumns(getsequence(small, 4)) == 6

for aln in (fasta, fasta_na, pfam, pfam_na)
  @test size(aln) == (4,110)
  @test length(aln) == 440
  @test ncolumns(aln) == 110
  @test nsequences(aln) == 4
end
