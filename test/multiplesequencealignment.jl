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
const pfam = readpfam("./data/PF09645_full.stockholm")

print("""

Parse Fasta
-----------
""")
const fasta = readfasta("./data/PF09645_full.fasta")

# Are the parsers getting the same result?
for field in msa_fields
  if field != :annotations
    @eval @test fasta.$field == pfam.$field
  end
end

