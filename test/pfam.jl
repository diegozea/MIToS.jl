# using Base.Test
# using MIToS.Pfam
# using MIToS.PDB
# using MIToS.Information
# using PairwiseListMatrices
using ROCAnalysis

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

print("""

Test Mapping PDB/Pfam
=====================
""")

let msa = read(joinpath(pwd(), "data", "PF09645_full.stockholm"), Stockholm, generatemapping=true, useidcoordinates=true),
    map = msacolumn2pdbresidue(msa, "F112_SSV1/3-112", "2VQC", "A", "PF09645", ascii(joinpath(pwd(), "data", "2vqc.xml.gz"))),
    res = residuesdict(read(joinpath(pwd(), "data", "2VQC.xml"), PDBML), "1", "A", "ATOM", "*")

  #     -45              20 pdb
  #.....QTLNSYKMAEIMYKILEK  msa seq
  #     123456789012345678  msa col
  #     345678901234567890  uniprot 3-20
  #    ****              *
  #12345678901234567890123  ColMap

  @test_throws KeyError map[5]  # insert
  @test map[6]  == ""           # missing
  @test map[7]  == "4"
  @test map[8]  == "5"
  @test map[23] == "20"

  #.....QTLNSYKMAEIMYKILEKKGELTLEDILAQFEISVPSAYNIQRALKAICERHPDECEVQYKNRKTTFKWIKQEQKEEQKQEQTQDNIAKIFDAQPANFEQTDQGFIKAKQ..... msa seq
  #.....X---HHHHHHHHHHHHHHHSEE-HHHHHHHH---HHHHHHHHHHHHHHHHH-TTTEEEEE-SS-EEEEE--XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX..... pdb ss
  #                                                                                                   111111111111111111111 ColMap hundreds
  #         111111111122222222223333333333444444444455555555556666666666777777777788888888889999999999000000000011111111112 ColMap tens
  #123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890 ColMap ones
  #                                                                           **                                     **

  @test_throws KeyError map[116]   # insert
  @test map[115] == ""             # missing
  @test map[77]  == ""             # missing
  @test map[76]  == "73"

  print("""

  Test Residues
  -------------
  """)

  msares = msaresidues(msa, res, map)

  @test_throws KeyError msares[5]   # insert
  @test_throws KeyError msares[6]   # missing
  @test msares[7].id.name  == "THR" # T
  @test msares[8].id.name  == "LEU" # L
  @test msares[23].id.name == "LYS" # K

  @test_throws KeyError msares[116]   # insert
  @test_throws KeyError msares[117]   # missing
  @test_throws KeyError msares[77]    # missing
  @test msares[76].id.name == "LYS"   # K

end

print("""

Test Contacts
=============
""")

let msa = read(joinpath(pwd(), "data", "PF09645_full.stockholm"), Stockholm, generatemapping=true, useidcoordinates=true),
    map = msacolumn2pdbresidue(msa, "F112_SSV1/3-112", "2VQC", "A", "PF09645", ascii(joinpath(pwd(), "data", "2vqc.xml.gz"))),
    res = residuesdict(read(joinpath(pwd(), "data", "2VQC.xml"), PDBML), "1", "A", "ATOM", "*")

  @test length(res) == 70
  @test length(unique(values(map))) == 71

  #     -45              20 pdb
  #.....QTLNSYKMAEIMYKILEK  msa seq
  #     123456789012345678  msa col
  #12345678901234567890123  ColMap
  #     345678901234567890  uniprot 3-20
  #     ***              *

  @test_throws KeyError res["3"]
  @test res[map[7]].id.number == "4"

  contacts = msacontacts(msa, res, map, 6.05)
  missings = sum(isnan(contacts), 1)

  @test size(contacts) == (110,110)

  @test missings[1]  == 110              # missing
  @test missings[2]  == (110 - (70 - 1)) # 70 residues, 1 diagonal NaN
  @test missings[3]  == (110 - (70 - 1))
  @test missings[18] == (110 - (70 - 1))

  @test missings[110] == 110                 # missing
  @test missings[72]  == 110                 # missing
  @test missings[71]  == (110 - (70 - 1))

  ncontacts = sum(contacts .== 1.0, 1)

  @test ncontacts[1]  == 0
  @test ncontacts[2]  == 2
  @test ncontacts[3]  == 6

  print("""

  Test MSA contact map using PDBResidues from the MSA
  ---------------------------------------------------
  """)

  println("")

  msares = msaresidues(msa, res, map)

  ncol = ncolumns(msa)
  colmap = getcolumnmapping(msa)

  for i in 1:(ncol-1), j in (i+1):ncol
    if (colmap[i] in keys(msares)) && (colmap[j] in keys(msares))
      print(".")
      @test (contacts[i,j] .== 1.0) == contact(msares[ colmap[i] ], msares[ colmap[j] ], 6.05)
    end
  end

  println("")

  print("""

  Test hasresidues
  """)

  mask = hasresidues(msa, map)

  @test mask[1] == false
  @test mask[2] == true
  @test sum(mask) == length(msares)

end

print("""

Test AUC and Contact Masks
==========================
""")

let ntru = 90,
    nfal = 100,
    score_tru = Float16[ 2 + 2x for x in randn(ntru)],
    score_fal = Float16[-2 + 2x for x in randn(nfal)],
    msacontacts = PairwiseListMatrix(vcat(ones(Float16, ntru), zeros(Float16, nfal)), Int16[x for x in 1:20], false),
    score = PairwiseListMatrix(vcat(score_tru, score_fal), Int16[x for x in 1:20], false),
    correct = 1 - auc(roc(score_tru, score_fal))

  @test AUC(score, msacontacts) == correct

end

print("""

AUC Example using MI
--------------------
""")

let msa = read(joinpath(pwd(), "data", "PF09645_full.stockholm"), Stockholm, generatemapping=true, useidcoordinates=true),
    map = msacolumn2pdbresidue(msa, "F112_SSV1/3-112", "2VQC", "A", "PF09645", ascii(joinpath(pwd(), "data", "2vqc.xml.gz"))),
    res = residuesdict(read(joinpath(pwd(), "data", "2VQC.xml"), PDBML), "1", "A", "ATOM", "*"),
    contacts = msacontacts(msa, res, map)

  @test round( AUC(buslje09(msa, lambda=0.05, threshold=62.0, samples=0)[2], contacts) , 4) == 0.5291
end
