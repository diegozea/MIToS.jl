using Base.Test
using MIToS.Pfam

print("""

Test download from Pfam
=======================
""")

let pfam_code = "PF11591"
  @test_throws ErrorException downloadpfam("2vqc")
  filename = downloadpfam(pfam_code)
  try
    aln = read(filename, Stockholm)
    if size(aln) == (6,34)
      @test getannotfile(aln, "ID") == "2Fe-2S_Ferredox"
    end
  finally
    rm(filename)
  end
end

print("""

Test PDB code from Pfam
=======================
""")

@test getseq2pdb(read(joinpath(pwd(), "data", "PF09645_full.stockholm"), Stockholm))["F112_SSV1/3-112"] == [("2VQC","A")]
