using Base.Test
using MIToS.Pfam
using MIToS.PDB

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
    map = msacolumn2pdbresidue("F112_SSV1/3-112", "2VQC", "A", "PF09645", msa, ascii(joinpath(pwd(), "data", "2vqc.xml.gz")))

  #     -45              20 pdb
  #.....QTLNSYKMAEIMYKILEK  msa seq
  #     123456789012345678  msa col
  #     345678901234567890  uniprot 3-20
  #    ****              *

  @test_throws KeyError map[0]  # insert
  @test map[1]  == ""           # missing
  @test map[2]  == "4"
  @test map[3]  == "5"
  @test map[18] == "20"

  #.....QTLNSYKMAEIMYKILEKKGELTLEDILAQFEISVPSAYNIQRALKAICERHPDECEVQYKNRKTTFKWIKQEQKEEQKQEQTQDNIAKIFDAQPANFEQTDQGFIKAKQ..... msa seq
  #.....X---HHHHHHHHHHHHHHHSEE-HHHHHHHH---HHHHHHHHHHHHHHHHH-TTTEEEEE-SS-EEEEE--XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX..... pdb ss
  #                                                                                                        11111111111      msa col hundreds
  #              11111111112222222222333333333344444444445555555555666666666677777777778888888888999999999900000000001      msa col tens
  #     12345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890      msa col ones
  #                                                                           **                                     **

  @test_throws KeyError map[111]   # insert
  @test map[110] == ""             # missing
  @test map[72]  == ""             # missing
  @test map[71]  == "73"

end

print("""

Test Contacts
=============
""")

let msa = read(joinpath(pwd(), "data", "PF09645_full.stockholm"), Stockholm, generatemapping=true, useidcoordinates=true),
    map = msacolumn2pdbresidue("F112_SSV1/3-112", "2VQC", "A", "PF09645", msa, ascii(joinpath(pwd(), "data", "2vqc.xml.gz"))),
    res = residuesdict(read(joinpath(pwd(), "data", "2VQC.xml"), PDBML), "1", "A", "ATOM", "*")

  @test length(res) == 70
  @test length(unique(values(map))) == 71

  #     -45              20 pdb
  #.....QTLNSYKMAEIMYKILEK  msa seq
  #     123456789012345678  msa col
  #     345678901234567890  uniprot 3-20
  #    ****              *

  @test_throws KeyError res["3"]
  @test res[map[2]].id.number == "4"

  contacts = msacontacts(res, map)
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

end
