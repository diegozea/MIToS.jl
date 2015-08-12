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
const F112_SSV1 = collect(".....QTLNSYKMAEIMYKILEKKGELTLEDILAQFEISVPSAYNIQRALKAICERHPDECEVQYKNRKTTFKWIKQEQKEEQKQEQTQDNIAKIFDAQPANFEQTDQGFIKAKQ.....")

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


