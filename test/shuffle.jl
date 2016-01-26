# using Base.Test
# using MIToS.MSA

print("""

Test shuffle_...
================
""")

let aln = read(joinpath(pwd(), "data", "PF09645_full.fasta.gz"), FASTA), copy_aln_msa = copy(aln.msa),
  gaps = (aln .== GAP), lcol = mean(aln .== Residue('L'), 1), lseq = mean(aln .== Residue('L'), 2)

  shuffle_residues_columnwise!(aln)
  @test (aln .== GAP) == gaps
  @test aln.msa != copy_aln_msa
  @test lcol == mean(aln .== Residue('L'), 1)
  @test lseq != mean(aln .== Residue('L'), 2)
end

let aln = read(joinpath(pwd(), "data", "PF09645_full.fasta.gz"), FASTA), copy_aln_msa = copy(aln.msa),
  gaps = (aln .== GAP), lcol = mean(aln .== Residue('L'), 1), lseq = mean(aln .== Residue('L'), 2)

  shuffle_residues_sequencewise!(aln)
  @test (aln .== GAP) == gaps
  @test aln.msa != copy_aln_msa
  @test lcol != mean(aln .== Residue('L'), 1)
  @test lseq == mean(aln .== Residue('L'), 2)
end

let aln = read(joinpath(pwd(), "data", "PF09645_full.fasta.gz"), FASTA), copy_aln_msa = copy(aln.msa),
  gaps = (aln .== GAP), lcol = mean(aln .== Residue('L'), 1), lseq = mean(aln .== Residue('L'), 2)

  shuffle_columnwise!(aln)
  @test (aln .== GAP) != gaps
  @test aln.msa != copy_aln_msa
  @test lcol == mean(aln .== Residue('L'), 1)
  @test lseq != mean(aln .== Residue('L'), 2)
end

let aln = read(joinpath(pwd(), "data", "PF09645_full.fasta.gz"), FASTA), copy_aln_msa = copy(aln.msa),
  gaps = (aln .== GAP), lcol = mean(aln .== Residue('L'), 1), lseq = mean(aln .== Residue('L'), 2)

  shuffle_sequencewise!(aln)
  @test (aln .== GAP) != gaps
  @test aln.msa != copy_aln_msa
  @test lcol != mean(aln .== Residue('L'), 1)
  @test lseq == mean(aln .== Residue('L'), 2)
end
