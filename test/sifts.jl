using Base.Test
using MIToS.SIFTS

print("""
Tests for SIFTS Mappings
========================
""")

print("""
Parse SIFTS
-----------
""")

const sifts_file = "./data/2vqc.xml.gz"

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
  @test_throws KeyError map[PDBeCoordinate(9)]  # Missing
  @test_throws KeyError map[PDBeCoordinate(80)] # Missing
  @test_throws KeyError map[PDBeCoordinate(1)]  # Missing
  @test map[PDBeCoordinate(10)] == PDBresnumCoordinate(4)
  @test map[PDBeCoordinate(79)] == PDBresnumCoordinate(73, "")
end

let map = siftsmapping(sifts_file, dbPDBe, "2vqc", dbPDB, "2vqc", chain="A")
  @test map[PDBeCoordinate(9)] == PDBresnumCoordinate(3)   # Missing
  @test map[PDBeCoordinate(80)] == PDBresnumCoordinate(74) # Missing
  @test map[PDBeCoordinate(1)] == PDBresnumCoordinate(-5)  # Missing # Negative Resnum
  @test map[PDBeCoordinate(10)] == PDBresnumCoordinate(4)
  @test map[PDBeCoordinate(79)] == PDBresnumCoordinate(73, "")
end

print("""
1SSX => Residues with insert codes: 15A 15B
""")

let map = siftsmapping("./data/1ssx.xml.gz", dbPDBe, "1ssx", dbPDB, "1ssx", chain="A")
  residue_A = map[PDBeCoordinate(1)]
  residue_B = map[PDBeCoordinate(2)]
  residue_C = map[PDBeCoordinate(3)]
  @test residue_A == PDBresnumCoordinate("15A")
  @test residue_B == PDBresnumCoordinate(15, 'B')
  @test residue_C == PDBresnumCoordinate(16)
  @test Int(residue_A) == 15
  @test Int(residue_B) == 15
  @test Int(residue_C) == 16
  @test residue_A.inscode == "A"
  @test residue_B.inscode == "B"
  @test residue_C.inscode == ""
end

let map = siftsmapping("./data/1ssx.xml.gz", dbPDB, "1ssx", dbUniProt, "P00778", chain="A")
  @test map[PDBresnumCoordinate("15A")] == UniProtCoordinate(200)
  @test map[PDBresnumCoordinate(15, 'B')] == UniProtCoordinate(201)
  @test map[PDBresnumCoordinate(16)] == UniProtCoordinate("202")
end

# print("""
# 1CBN => Multiple InterProt annotations for Thr2
# Since InterPro ResNum can be "None" instead of a Int; InterPro is not supported
# """)

# let map = siftsmapping("./data/1cbn.xml.gz", "A", dbPDBe, dbInterPro, "IPR001010")
#   @test_throws KeyError map[PDBeCoordinate(1)] # Without InterPro
#   @test map[PDBeCoordinate(2)] == UniProtCoordinate(2) # Same ResNum for different InterPros
# end

# let map = siftsresidues("./data/1cbn.xml.gz", chain="A")
#   @test length(map[1].InterPro) == 0 # Without InterPro
#   @test length(map[2].InterPro) == 4 # Same ResNum for different InterPros
# end

print("""
1AS5 => NMR
""")

let map = siftsmapping("./data/1as5.xml.gz", dbPDBe,"1as5", dbUniProt, "P56529", chain="A") # missings=true : NMR there are not missing residues
  @test map[PDBeCoordinate(23)] == UniProtCoordinate(73)
  @test_throws KeyError map[PDBeCoordinate(24)] # Without UniProt
end

let map = siftsmapping("./data/1as5.xml.gz", dbPDBe,"1as5", dbPDB, "1as5", chain="A")
  @test map[PDBeCoordinate(24)] == PDBresnumCoordinate(24)
end

print("""
1DPO => Inserted residues lack insertion letters
Single unnamed chain in 1DPO contains insertions at postions 184 (Gly, Phe), 188 (Gly, Lys), and 221 (Ala, Leu) but no insertion letters.
""")

let map = siftsmapping("./data/1dpo.xml.gz", dbPDBe, "1dpo", dbPDB, "1dpo", chain="A") # Unnamed chain is "A" in SIFTS
  @test map[PDBeCoordinate(164)] == PDBresnumCoordinate(184)
  @test map[PDBeCoordinate(165)] == PDBresnumCoordinate(184, 'A') # Has insertion code in SIFTS
  @test map[PDBeCoordinate(169)] == PDBresnumCoordinate(188)
  @test map[PDBeCoordinate(170)] == PDBresnumCoordinate(188, 'A') # Has insertion code in SIFTS
  @test map[PDBeCoordinate(198)] == PDBresnumCoordinate(221)
  @test map[PDBeCoordinate(199)] == PDBresnumCoordinate(221, 'A') # Has insertion code in SIFTS
end

print("""
1IGY => Insertions have more than one copy of the same amino acid in a single insertion block.
For example, chain B in 1IGY contains a block of four residues inserted at sequence position 82. The block contains Leu-Ser-Ser-Leu.
""")

let map = siftsmapping("./data/1igy.xml.gz", dbPDBe, "1igy", dbCATH, "2.60.40.10", chain="B")
  @test map[PDBeCoordinate(82)] == PDBresnumCoordinate(82)
  @test map[PDBeCoordinate(83)] == PDBresnumCoordinate("82A")
  @test map[PDBeCoordinate(84)] == PDBresnumCoordinate("82B")
  @test map[PDBeCoordinate(85)] == PDBresnumCoordinate("82C")
end

print("""
1HAG => Chain E begins with 1H, 1G, 1F, ... 1A, then 1 (in reverse alphabetic order)
""")

let map = siftsmapping("./data/1hag.xml.gz", dbPDBe, "1hag", dbPDB, "1hag", chain="E")
  @test map[PDBeCoordinate(1)] == PDBresnumCoordinate(1, 'H')
  @test map[PDBeCoordinate(2)] == PDBresnumCoordinate(1, 'G')
  @test map[PDBeCoordinate(3)] == PDBresnumCoordinate(1, 'F')
  @test map[PDBeCoordinate(4)] == PDBresnumCoordinate(1, 'E')
  @test map[PDBeCoordinate(5)] == PDBresnumCoordinate(1, 'D')
  @test map[PDBeCoordinate(6)] == PDBresnumCoordinate(1, 'C')
  @test map[PDBeCoordinate(7)] == PDBresnumCoordinate(1, 'B')
  @test map[PDBeCoordinate(8)] == PDBresnumCoordinate(1, 'A')
  @test map[PDBeCoordinate(9)] == PDBresnumCoordinate(1)
end

print("""

Functions for SIFTSResidues
===========================
""")

print("""
1NSA => Contains a single (unnamed) protein chain with sequence 7A-95A that continues 4-308.
""")

let map = siftsresidues("./data/1nsa.xml.gz")
  four = findfirst(x -> has(x, dbPDB, "1nsa", PDBresnumCoordinate("4")), map)
  @test findfirst(x -> has(x, dbPDB, "1nsa", PDBresnumCoordinate("95A")), map) + 1 == four
  @test getcoordinate(map[four], dbPDB, "1nsa") == PDBresnumCoordinate("4")
end

print("""
1IAO => Contains in chain B (in this order) 1S, 323P-334P, 6-94, 94A, 95-188, 1T, 2T
""")

let map = siftsresidues("./data/1iao.xml.gz")
  i = findfirst(x -> has(x, dbPDB, "1iao", PDBresnumCoordinate("1S")) && ischain(x, "B"), map)
  getcoordinate(map[i+2], dbPDB, "1iao", "B") == PDBresnumCoordinate("323P")
end
