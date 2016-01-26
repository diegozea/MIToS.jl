# using Base.Test
# using MIToS.MSA

# Fields of MultipleSequenceAlignment
const msa_fields = Symbol[:id, :msa, :annotations]

# Fields of AlignedAlignedSequence
const seq_fields = Symbol[:id, :index, :sequence, :annotations]

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

const pfam = read(joinpath(pwd(), "data", "PF09645_full.stockholm"), Stockholm)
const pfam_na = read(joinpath(pwd(), "data", "PF09645_full.stockholm"), Stockholm, MultipleSequenceAlignment)
const F112_SSV1 = collect(".....QTLNSYKMAEIMYKILEKKGELTLEDILAQFEISVPSAYNIQRALKAICERHPDECEVQYKNRKTTFKWIKQEQKEEQKQEQTQDNIAKIFDAQPANFEQTDQGFIKAKQ.....")

let rawpfam = read(joinpath(pwd(), "data", "PF09645_full.stockholm"), Stockholm, Matrix{Residue})
  @test rawpfam == pfam.msa
end

@test pfam_na == convert(MultipleSequenceAlignment, pfam)

@test size(pfam.msa, 1) == 4
@test size(pfam.msa, 2)  == length(F112_SSV1[ F112_SSV1 .!= '.' ]) # Without inserts
@test slice(pfam.msa,4,:) == convert(Vector{Residue}, F112_SSV1[ F112_SSV1 .!= '.' ])
@test pfam.id.items == ["C3N734_SULIY/1-95", "H2C869_9CREN/7-104", "Y070_ATV/2-70", "F112_SSV1/3-112"]
@test !isempty(pfam.annotations)
@test length(pfam.annotations.columns) == 2
@test length(pfam.annotations.residues) == 1
@test pfam.annotations.residues[("F112_SSV1/3-112","SS")] == "X---HHHHHHHHHHHHHHHSEE-HHHHHHHH---HHHHHHHHHHHHHHHHH-TTTEEEEE-SS-EEEEE--XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
@test pfam.annotations.columns["seq_cons"] == "...NshphAclhaKILppKtElolEDIlAQFEISsosAYsI.+sL+hICEpH.-ECpsppKsRKTlhh.hKpEphppptpEp..ppItKIhsAp................"
@test pfam.annotations.sequences[("F112_SSV1/3-112","DR")] == "PDB; 2VQC A; 4-73;"

print("""
Test parse for string inputs
""")

const pfam_string = """# STOCKHOLM 1.0
#=GS C3N734_SULIY/1-95   AC C3N734.1
#=GS H2C869_9CREN/7-104  AC H2C869.1
#=GS Y070_ATV/2-70       AC Q3V4T1.1
#=GS F112_SSV1/3-112     AC P20220.1
#=GS F112_SSV1/3-112     DR PDB; 2VQC A; 4-73;
C3N734_SULIY/1-95              ...mp---NSYQMAEIMYKILQQKKEISLEDILAQFEISASTAYNVQRTLRMICEKHPDECEVQTKNRRTIFKWIKNEETTEEGQEE--QEIEKILNAQPAE-------------k....
H2C869_9CREN/7-104             ...nk--LNDVQRAKLLVKILQAKGELDVYDIMLQFEISYTRAIPIMKLTRKICEAQ-EICTYDEKEHKLVSLNAKKEKVEQDEEENEREEIEKILDAH----------------trreq
Y070_ATV/2-70                  qsvne-------VAQQLFSKLREKKEITAEDIIAIYNVTPSVAYAIFTVLKVMCQQHQGECQAIKRGRKTVI-------------------------------------------vskq.
F112_SSV1/3-112                .....QTLNSYKMAEIMYKILEKKGELTLEDILAQFEISVPSAYNIQRALKAICERHPDECEVQYKNRKTTFKWIKQEQKEEQKQEQTQDNIAKIFDAQPANFEQTDQGFIKAKQ.....
#=GR F112_SSV1/3-112     SS    .....X---HHHHHHHHHHHHHHHSEE-HHHHHHHH---HHHHHHHHHHHHHHHHH-TTTEEEEE-SS-EEEEE--XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX.....
#=GC SS_cons                   .....X---HHHHHHHHHHHHHHHSEE-HHHHHHHH---HHHHHHHHHHHHHHHHH-TTTEEEEE-SS-EEEEE--XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX.....
#=GC seq_cons                  ........NshphAclhaKILppKtElolEDIlAQFEISsosAYsI.+sL+hICEpH.-ECpsppKsRKTlhh.hKpEphppptpEp..ppItKIhsAp................h....
//
"""

@test parse(pfam_string, Stockholm) == pfam

print("""

Parse Fasta
-----------
""")

const fasta = read(joinpath(pwd(), "data", "PF09645_full.fasta.gz"), FASTA)
const small = read(joinpath(pwd(), "data", "Gaoetal2011.fasta"), FASTA)
const fasta_na = read(joinpath(pwd(), "data", "PF09645_full.fasta.gz"), FASTA, MultipleSequenceAlignment)
const small_na = read(joinpath(pwd(), "data", "Gaoetal2011.fasta"), FASTA, MultipleSequenceAlignment, deletefullgaps=false)

let rawfasta = read(joinpath(pwd(), "data", "PF09645_full.fasta.gz"), FASTA, Matrix{Residue})
  @test rawfasta == fasta.msa
end

@test fasta_na == convert(MultipleSequenceAlignment, fasta)
@test small_na == convert(MultipleSequenceAlignment, small)

print("""
Test download of read(URL, ...)
""")

@test read("https://raw.githubusercontent.com/diegozea/MIToS.jl/master/test/data/Gaoetal2011.fasta", FASTA) == small
@test read("https://raw.githubusercontent.com/diegozea/MIToS.jl/master/test/data/PF09645_full.fasta.gz", FASTA) == fasta

print("""
Test parse for string inputs
""")

const fasta_string = """>SEQ1
DAWAEE
>SEQ2
DAWAEF
>SEQ3
DAWAED
>SEQ4
DAYCMD
>SEQ5
DAYCMT
>SEQ6
DAYCMT
"""

@test parse(fasta_string, FASTA) == small

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
@test small.id.items == ASCIIString[ "SEQ$i" for i in 1:6 ]

print("""
Are the parsers (pfam/fasta) getting the same result?
""")

for field in msa_fields
  if field != :annotations
    @eval @test fasta.$field == pfam.$field
  end
end

print("""
Test parse with ambiguous or not standard residues (are gaps on MIToS)
""")

let default = read(joinpath(pwd(), "data", "alphabet.fasta"), FASTA, generatemapping=true),
    notused = read(joinpath(pwd(), "data", "alphabet.fasta"), FASTA, checkalphabet=true) # deletenotalphabetsequences! is called also for Raw & Stockholm

  @test Base.vec(default[1,:]) == res"ARNDCQEGHILKMFPSTWYV"
  for i in 2:nsequences(default)
    @test Base.vec(default[i,:]) == res"AR-DCQEGHILKMFPSTWYV"
    @test  getsequencemapping(default,1) == getsequencemapping(default,i)
  end

  @test nsequences(notused) == 1
  @test nsequences(read(joinpath(pwd(), "data", "alphabet.fasta"), FASTA, Matrix{Residue}, checkalphabet=true)) == 1 # _strings_to_msa is called also for Raw & Stockholm
end

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
const copy_small = copy(small)
const deepcopy_small_na = deepcopy(small_na)

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

print("""

AnnotatedAlignedSequence and AlignedSequence
--------------------------------------------
""")

print("""
Test getsequence using Index and ID and convert to AlignedSequence
""")

for i in 1:4
  annseq = getsequence(pfam, i)
  @test isa(annseq, AnnotatedAlignedSequence)
  seq = getsequence(pfam_na, i)
  @test isa(seq, AlignedSequence)
  @test AlignedSequence(annseq) == seq
end

for id in ["C3N734_SULIY/1-95", "H2C869_9CREN/7-104", "Y070_ATV/2-70", "F112_SSV1/3-112"]
  annseq = getsequence(pfam, id)
  @test pfam[id] == annseq
  @test isa(annseq, AnnotatedAlignedSequence)
  seq = getsequence(pfam_na, id)
  @test pfam_na[id] == seq
  @test isa(seq, AlignedSequence)
  @test AlignedSequence(annseq) == seq
end

print("""

Test Filters!
-------------
""")

print("""
filtersequences! with and without annotations
""")
@test getsequence(filtersequences!(copy(pfam), [1,2,3,4] .== 3), 1) == getsequence(pfam,3)
@test getsequence(filtersequences!(copy(pfam), [1,2,3,4] .> 2), 2) == getsequence(pfam,4)

@test getsequence(filtersequences!(copy(pfam_na), [1,2,3,4] .== 3), 1) == getsequence(pfam_na,3)
@test getsequence(filtersequences!(copy(pfam_na), [1,2,3,4] .> 2), 2) == getsequence(pfam_na,4)

print("""
filtercolumns! with and without annotations
""")

@test getresidues(getsequence(filtercolumns!(copy(pfam), collect(1:110) .<= 10), 4)) == res"QTLNSYKMAE"
@test getresidues(getsequence(filtercolumns!(copy(pfam_na), collect(1:110) .<= 10), 4)) == res"QTLNSYKMAE"

print("""
Test for BoundsError's on filters with bad masks
""")

@test_throws BoundsError filtersequences!(copy(pfam), [1,2,3,4,5,6,7,8,9,10] .> 2)
@test_throws BoundsError filtersequences!(copy(pfam_na), [1,2,3] .> 2)
@test_throws BoundsError filtercolumns!(copy(pfam), [1,2,3] .> 2)
@test_throws BoundsError filtercolumns!(copy(pfam_na), collect(1:200) .<= 10)

print("""
Test copy and deepcopy for sequences
""")

for (index, value) in [ (2,Residue('H')) , (2:3,res"HG") ]
  annseq = getsequence(pfam, 4)
  seq = getsequence(pfam_na, 4)
  copied_annseq = copy(annseq)
  deepcopied_seq = deepcopy(seq)
  copied_annseq[index] = value
  deepcopied_seq[index] = value
  @test copied_annseq[index] == value
  @test annseq[index] != value
  @test deepcopied_seq[index] == value
  @test seq[index] != value
end

print("""
Test filtercolumns! for sequences
""")

let i = 4
  annseq = getsequence(pfam, i)
  seq = getsequence(pfam_na, i)
  print("""
  Test Sequence .== Residue for mask generation
  """)
  @test annseq[annseq .== Residue('Q')] == res"QQQQQQQQQQQQQQ"
  @test seq[seq .== Residue('Q')] == res"QQQQQQQQQQQQQQ"
  print("""
  Test filtercolumns!
  """)
  filtered_annseq = filtercolumns!(copy(annseq), annseq .== Residue('Q'))
  filtered_seq = filtercolumns!(copy(seq), seq .== Residue('Q'))
  @test_throws BoundsError filtercolumns!(copy(annseq), collect(1:(length(seq)-10)) .> 2)
  @test_throws BoundsError filtercolumns!(copy(seq), collect(1:(length(seq)+10)) .<= 10)
  @test getresidues(filtered_annseq) == res"QQQQQQQQQQQQQQ"
  @test getresidues(filtered_seq) == res"QQQQQQQQQQQQQQ"
  @test AlignedSequence(filtered_annseq) == filtered_seq
  @test getannotcolumn(filtered_annseq, "SS_cons") == "XHHEXXXXXXXXXX"
  @test getannotresidue(filtered_annseq, "F112_SSV1/3-112", "SS") == "XHHEXXXXXXXXXX"
end

print("""

Test Annotation Getters
-----------------------
""")

@test getannotcolumn(pfam, "SS_cons") == getannotcolumn(getsequence(pfam, 4), "SS_cons")
@test getannotresidue(pfam, "F112_SSV1/3-112", "SS") == getannotresidue(getsequence(pfam, 4), "F112_SSV1/3-112", "SS")
@test getannotfile(pfam) == getannotfile(getsequence(pfam, 4))
@test getannotsequence(pfam,"F112_SSV1/3-112","DR") == getannotsequence(getsequence(pfam, 4), "F112_SSV1/3-112", "DR")

print("""

Test Gaps and Coverage
----------------------
""")

@test gappercentage(small) == 0.0
@test gappercentage(getsequence(small,1)) == 0.0
@test gappercentage(small[1,:]) == 0.0

@test gappercentage(getsequence(pfam,4)) == 0.0
@test gappercentage(getsequence(pfam_na,4)) == 0.0

@test gappercentage(convert(Vector{Residue}, F112_SSV1)) == mean(F112_SSV1 .== '.')

@test gappercentage(pfam[:,1]) == 0.75
@test gappercentage(pfam_na[:,1]) == 0.75

print("""
Test residuepercentage
""")

@test residuepercentage(small) == 1.0
@test residuepercentage(getsequence(pfam, 4)) == 1.0
@test residuepercentage(pfam[:,1]) == 0.25

print("""
Test coverage
""")

@test coverage(pfam)[4] == 1.0
@test sum(coverage(small)) == 6.0

print("""

Reference
---------
""")

print("""
Test setreference!
""")

let copy_pfam = copy(pfam), copy_pfam_na = copy(pfam_na)
  setreference!(copy_pfam, 4)
  @test getresidues(getsequence(copy_pfam,1)) == res"QTLNSYKMAEIMYKILEKKGELTLEDILAQFEISVPSAYNIQRALKAICERHPDECEVQYKNRKTTFKWIKQEQKEEQKQEQTQDNIAKIFDAQPANFEQTDQGFIKAKQ"
  setreference!(copy_pfam, "C3N734_SULIY/1-95")
  @test getresidues(getsequence(copy_pfam,4)) == res"QTLNSYKMAEIMYKILEKKGELTLEDILAQFEISVPSAYNIQRALKAICERHPDECEVQYKNRKTTFKWIKQEQKEEQKQEQTQDNIAKIFDAQPANFEQTDQGFIKAKQ"
  @test copy_pfam == pfam

  setreference!(copy_pfam_na, 4)
  @test getresidues(getsequence(copy_pfam_na,1)) == res"QTLNSYKMAEIMYKILEKKGELTLEDILAQFEISVPSAYNIQRALKAICERHPDECEVQYKNRKTTFKWIKQEQKEEQKQEQTQDNIAKIFDAQPANFEQTDQGFIKAKQ"
  setreference!(copy_pfam_na, "C3N734_SULIY/1-95")
  @test getresidues(getsequence(copy_pfam_na,4)) == res"QTLNSYKMAEIMYKILEKKGELTLEDILAQFEISVPSAYNIQRALKAICERHPDECEVQYKNRKTTFKWIKQEQKEEQKQEQTQDNIAKIFDAQPANFEQTDQGFIKAKQ"
  @test copy_pfam_na == pfam_na
end

print("""
Test gapstrip!
""")

let copy_pfam = copy(pfam), copy_pfam_na = copy(pfam_na)
  setreference!(copy_pfam, 4)
  setreference!(copy_pfam_na, 4)
  gapstrip!(copy_pfam, gaplimit=1.0, coveragelimit=0.0)
  gapstrip!(copy_pfam_na, gaplimit=1.0, coveragelimit=0.0)
  @test size(copy_pfam) == (4, 110)
  @test copy_pfam == copy_pfam_na

  setreference!(copy_pfam, 1)
  setreference!(copy_pfam_na, 1)
  gapstrip!(copy_pfam)
  gapstrip!(copy_pfam_na)
  residuepercentage(copy_pfam[1,:]) == 1.0
  residuepercentage(copy_pfam_na[1,:]) == 1.0
end

print("""
Test adjustreference!
""")

let copy_pfam = copy(pfam), copy_pfam_na = copy(pfam_na)
  adjustreference!(copy_pfam)
  adjustreference!(copy_pfam_na)
  residuepercentage(copy_pfam[1,:]) == 1.0
  residuepercentage(copy_pfam_na[1,:]) == 1.0
end

print("""

Printers
--------
""")

print("""
Test asciisequence
""")

@test asciisequence(small, 4) == "DAYCMD"

print("""
Test printfasta
""")

let io = IOBuffer()
  print(io, small, FASTA)
  @test takebuf_string(io) == fasta_string
end

print("""
Test printpfam
""")

let io = IOBuffer()
  print(io, pfam, Stockholm)
  @test parse(takebuf_string(io), Stockholm) == pfam
end

let io = IOBuffer()
  print(io, small_na, Stockholm)
  @test parse(takebuf_string(io), Stockholm) == small_na
end
