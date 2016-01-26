# using Base.Test
# using MIToS.SIFTS

print("""

Tests for SIFTS Mappings
========================
""")

print("""

Parse SIFTS
-----------
""")

const sifts_file = joinpath(pwd(), "data", "2vqc.xml.gz")

print("""
all against all
""")

let dbs = [(dbUniProt, "P20220"), (dbPfam, "PF09645"), (dbNCBI, "244589"), (dbPDB, "2vqc"), (dbSCOP, "153426"), (dbPDBe, "2vqc")]
  for (to,toid) in dbs, (from, fromid) in dbs
    map = siftsmapping(sifts_file, from, fromid, to, toid)
    for (k,v) in map
      to == from && @test k == v
    end
  end
end

print("""

siftsmapping
------------
""")

print("""
missings
""")

let map = siftsmapping(sifts_file, dbPDBe, "2vqc", dbPDB, "2vqc", chain="A", missings=false)
  @test_throws KeyError map[9]  # Missing
  @test_throws KeyError map[80] # Missing
  @test_throws KeyError map[1]  # Missing
  @test map[10] == "4"
  @test map[79] == "73"
end

let map = siftsmapping(sifts_file, dbPDBe, "2vqc", dbPDB, "2vqc", chain="A")
  @test map[9] == "3"   # Missing
  @test map[80] == "74" # Missing
  @test map[1] == "-5"  # Missing # Negative Resnum
  @test map[10] == "4"
  @test map[79] == "73"
end

print("""
1SSX => Residues with insert codes: 15A 15B
""")

let map = siftsmapping(joinpath(pwd(), "data", "1ssx.xml.gz"), dbPDBe, "1ssx", dbPDB, "1ssx", chain="A")
  residue_A = map[1]
  residue_B = map[2]
  residue_C = map[3]
  @test residue_A == "15A"
  @test residue_B == "15B"
  @test residue_C == "16"
end

let map = siftsmapping(joinpath(pwd(), "data", "1ssx.xml.gz"), dbPDB, "1ssx", dbUniProt, "P00778", chain="A")
  @test map["15A"] == 200
  @test map["15B"] == 201
  @test map["16"] == 202
end

print("""
1CBN => Multiple InterProt annotations, the last is used.
Identical PDBe ResNum for Residue 22:

        <residue dbSource="PDBe" dbCoordSys="PDBe" dbResNum="22" dbResName="SER">
          <crossRefDb dbSource="PDB" dbCoordSys="PDBresnum" dbAccessionId="1cbn" dbResNum="22" dbResName="SER" dbChainId="A"/>
          <crossRefDb dbSource="UniProt" dbCoordSys="UniProt" dbAccessionId="P01542" dbResNum="22" dbResName="P"/>
          ....
        <residue dbSource="PDBe" dbCoordSys="PDBe" dbResNum="22" dbResName="PRO">
          <crossRefDb dbSource="PDB" dbCoordSys="PDBresnum" dbAccessionId="1cbn" dbResNum="22" dbResName="PRO" dbChainId="A"/>
          <crossRefDb dbSource="UniProt" dbCoordSys="UniProt" dbAccessionId="P01542" dbResNum="22" dbResName="P"/>
          ...

""")

let map = siftsmapping(joinpath(pwd(), "data", "1cbn.xml.gz"), dbPDBe, "1cbn", dbInterPro, "IPR001010", chain="A")
  @test_throws KeyError map[1] # Without InterPro
  @test map[2] == "2" # Same ResNum for different InterPros
end

let map = read(joinpath(pwd(), "data", "1cbn.xml.gz"), SIFTSXML, chain="A")
  @test length(map[1].InterPro) == 0 # Without InterPro
  @test length(map[2].InterPro) == 4 # Same ResNum for different InterPros
end

print("""
1AS5 => NMR
""")

let map = siftsmapping(joinpath(pwd(), "data", "1as5.xml.gz"), dbPDBe,"1as5", dbUniProt, "P56529", chain="A") # missings=true : NMR there are not missing residues
  @test map[23] == 73
  @test_throws KeyError map[24] # Without UniProt
end

let map = siftsmapping(joinpath(pwd(), "data", "1as5.xml.gz"), dbPDBe,"1as5", dbPDB, "1as5", chain="A")
  @test map[24] == "24"
end

print("""
1DPO => Inserted residues lack insertion letters
Single unnamed chain in 1DPO contains insertions at postions 184 (Gly, Phe), 188 (Gly, Lys), and 221 (Ala, Leu) but no insertion letters.
""")

let map = siftsmapping(joinpath(pwd(), "data", "1dpo.xml.gz"), dbPDBe, "1dpo", dbPDB, "1dpo", chain="A") # Unnamed chain is "A" in SIFTS
  @test map[164] == "184"
  @test map[165] == "184A" # Has insertion code in SIFTS
  @test map[169] == "188"
  @test map[170] == "188A" # Has insertion code in SIFTS
  @test map[198] == "221"
  @test map[199] == "221A" # Has insertion code in SIFTS
end

print("""
1IGY => Insertions have more than one copy of the same amino acid in a single insertion block.
For example, chain B in 1IGY contains a block of four residues inserted at sequence position 82. The block contains Leu-Ser-Ser-Leu.
""")

let map = siftsmapping(joinpath(pwd(), "data", "1igy.xml.gz"), dbPDBe, "1igy", dbCATH, "2.60.40.10", chain="B")
  @test map[82] == "82"
  @test map[83] == "82A"
  @test map[84] == "82B"
  @test map[85] == "82C"
end

print("""
1HAG => Chain E begins with 1H, 1G, 1F, ... 1A, then 1 (in reverse alphabetic order)
""")

let map = siftsmapping(joinpath(pwd(), "data", "1hag.xml.gz"), dbPDBe, "1hag", dbPDB, "1hag", chain="E")
  @test map[1] == "1H"
  @test map[2] == "1G"
  @test map[3] == "1F"
  @test map[4] == "1E"
  @test map[5] == "1D"
  @test map[6] == "1C"
  @test map[7] == "1B"
  @test map[8] == "1A"
  @test map[9] == "1"
end

print("""

Functions for SIFTSResidues
===========================
""")

print("""
1NSA => Contains a single (unnamed) protein chain with sequence 7A-95A that continues 4-308.
""")

print("""
Test findobjects & read
""")

let mapp = read(joinpath(pwd(), "data", "1nsa.xml.gz"), SIFTSXML)
  four = findobjects(mapp, dbPDB, Is(:id, "1nsa"), Is(:number, "4"))[1]
  @test findobjects(mapp, dbPDB, Is(:id, "1nsa"), Is(:number, "95A"))[1] + 1 == four
  @test get(mapp[ findobjects(mapp, dbPDB, Is(:id, "1nsa"), Is(:number, "95A")) ][1].PDB).number == "95A"
end

print("""
Test collectcaptures
""")

let cap = collectcaptures(read(joinpath(pwd(), "data", "1nsa.xml.gz"), SIFTSXML), dbPDB, :number, dbPDB, Is(:id, "1nsa"), Is(:number, "95A"))
  @test get(cap[ [!isnull(x) for x in cap] ][1]) == "95A"
end

print("""
1IAO => Contains in chain B (in this order) 1S, 323P-334P, 6-94, 94A, 95-188, 1T, 2T
""")

print("""
Test collectobjects and isobject
""")

let mapp = read(joinpath(pwd(), "data", "1iao.xml.gz"), SIFTSXML)
  @test collectobjects(mapp, dbPDB, Is(:id, "1iao"), Is(:number, "1S"), Is(:chain, "B"))[1].PDBe.number == 1
  i = findobjects(mapp, dbPDB, Is(:id, "1iao"), Is(:number, "1S"), Is(:chain, "B"))[1]
  @test isobject(mapp[i+2], dbPDB, Is(:id, "1iao"), Is(:chain, "B"), Is(:number, "323P"))
end

print("""

Test download
=============
""")

let pdb = "2vqc"
  filename = downloadsifts(pdb)
  try
    @test_throws ErrorException downloadsifts("2vqc_A")
    @test length(read(filename, SIFTSXML)) == length(read(joinpath(pwd(), "data", "$(pdb).xml.gz"), SIFTSXML))
  finally
    rm(filename)
  end
end
