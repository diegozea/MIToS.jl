using Base.Test
using MIToS.MSA

print("""

Tests for Multiple Sequence Alignments
======================================
""")

print("""

Parse Pfam
----------
""")
const pfam = readpfam("../test/data/PF09645_full.stockholm")

print("""

Parse Fasta
-----------
""")
const fasta = readfasta("../test/data/PF09645_full.fasta")

@test fasta.msa == pfam.msa
