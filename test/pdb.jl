# using Base.Test
# using MIToS.PDB

txt(code) = joinpath(pwd(), "data", string(uppercase(code), ".pdb"))
xml(code) = joinpath(pwd(), "data", string(uppercase(code), ".xml"))

print("""

Tests for PDB
=============
""")

print("""

Parse PDB and PDBML
-------------------
""")

print("""
2VQC => Missings
""")

let code = "2VQC"
  pdb = read(txt(code), PDBFile)
  pdbml = read(xml(code), PDBML)
  @test findfirst(x -> x.id.number == "4",  pdb) == findfirst(x -> x.id.number == "4",  pdbml)
  @test findfirst(x -> x.id.number == "73", pdb) == findfirst(x -> x.id.number == "73", pdbml)

end

print("""
Test download
""")

let code = "2VQC"
  pdb = read(txt(code), PDBFile)
  filename = downloadpdb(code)
  try
    pdbml = read(filename, PDBML)
    @test findfirst(x -> x.id.number == "4",  pdb) == findfirst(x -> x.id.number == "4",  pdbml)
    @test findfirst(x -> x.id.number == "73", pdb) == findfirst(x -> x.id.number == "73", pdbml)
  finally
    rm(filename)
  end
end

print("""
1SSX => Residues with insert codes: 15A 15B
""")

let code = "1SSX"
  pdb = read(txt(code), PDBFile)
  pdbml = read(xml(code), PDBML)

  @test findobjects(pdbml, Is(:number, "15A"))[1] == 1
  @test findobjects(pdbml, Is(:number, "15B"))[1] == 2

  @test findobjects(pdb, Is(:number, "15A"))[1] == 1
  @test findobjects(pdb, Is(:number, "15B"))[1] == 2

  print("""
  test @residues
  """)
  @test (@residues pdb model "*" chain "*" group "*" residue "141")[1] == collectobjects(pdb, Is(:number, "141"))[1]
  # Testing the macro in let block:
  mo = "1"
  ch = "A"
  gr = "*"
  re = "141"
  @test (@residues pdb model mo chain ch group gr residue re)[1] == (@residues pdb model "*" chain "*" group gr residue "141")[1]

  print("""
  Occupancy != 1.0 and @atom
  """)
  @test sum([ get(occ,0) for occ in  collectcaptures(collectobjects(pdbml, Is(:number, "141"))[1].atoms, :occupancy, Is(:atom, "HH22")) ]) == 1.0
  @test sum( [ atom.occupancy for atom in @atoms pdbml model "1" chain "A"  group "*" residue "141" atom "HH22" ] ) == 1.0
  # Testing the macro in let block:
  at = "HH22"
  @test sum( [ atom.occupancy for atom in @atoms pdbml model mo chain ch group gr residue re atom at ] ) == 1.0

end

print("""
1CBN => Identical PDBe ResNum for Residue 22:

        <residue dbSource="PDBe" dbCoordSys="PDBe" dbResNum="22" dbResName="SER">
          <crossRefDb dbSource="PDB" dbCoordSys="PDBresnum" dbAccessionId="1cbn" dbResNum="22" dbResName="SER" dbChainId="A"/>
          <crossRefDb dbSource="UniProt" dbCoordSys="UniProt" dbAccessionId="P01542" dbResNum="22" dbResName="P"/>
          ....
        <residue dbSource="PDBe" dbCoordSys="PDBe" dbResNum="22" dbResName="PRO">
          <crossRefDb dbSource="PDB" dbCoordSys="PDBresnum" dbAccessionId="1cbn" dbResNum="22" dbResName="PRO" dbChainId="A"/>
          <crossRefDb dbSource="UniProt" dbCoordSys="UniProt" dbAccessionId="P01542" dbResNum="22" dbResName="P"/>
          ...

""")

let code = "1CBN"
  pdb = read(txt(code), PDBFile)
  pdbml = read(xml(code), PDBML)

  @test length( @residues pdb   model "1" chain "A" group "*" residue "22" ) == 2
  @test length( @residues pdbml model "1" chain "A" group "*" residue "22" ) == 2

  @test ASCIIString[ res.id.name for res in  @residues pdb   model "1" chain "A" group "*" residue "22" ] == ["SER", "PRO"]
  @test ASCIIString[ res.id.name for res in  @residues pdbml model "1" chain "A" group "*" residue "22" ] == ["SER", "PRO"]
end

print("""
1AS5 => NMR
""")

let code = "1AS5"
  pdb = read(txt(code), PDBFile)
  pdbml = read(xml(code), PDBML)

  @test length( @residues pdbml model "1"  chain "A" group "*" residue "*" ) == 25
  @test length( @residues pdbml model "14" chain "A" group "*" residue "*" ) == 25

  @test length( @residues pdbml model "*"  chain "A" group "*" residue "*" ) == 25*14
end

print("""
1DPO => Inserted residues lack insertion letters
Single unnamed chain in 1DPO contains insertions at postions 184 (Gly, Phe), 188 (Gly, Lys), and 221 (Ala, Leu) but no insertion letters.
""")

let code = "1DPO"
  pdb = read(txt(code), PDBFile)
  pdbml = read(xml(code), PDBML)

  # Single "A" chain for PDB (auth_asym_id in PDBML)
  @test unique(ASCIIString[ res.id.chain for res in pdb  ]) == ASCIIString[ "A" ]
  # But 'A':'H' chains for PDBML (label_asym_id)
  @test unique(ASCIIString[ res.id.chain for res in pdbml]) == ASCIIString[ string(chain) for chain in 'A':'H' ]

  @test ASCIIString[ res.id.name for res in @residues pdb   model "1" chain "A" group "*" residue r"^184[A-Z]?$" ] == ASCIIString["GLY", "PHE"]
  @test ASCIIString[ res.id.name for res in @residues pdbml model "1" chain "A" group "*" residue r"^184[A-Z]?$" ] == ASCIIString["GLY", "PHE"]

  @test ASCIIString[ res.id.name for res in @residues pdb   model "1" chain "A" group "*" residue r"^188[A-Z]?$" ] == ASCIIString["GLY", "LYS"]
  @test ASCIIString[ res.id.name for res in @residues pdbml model "1" chain "A" group "*" residue r"^188[A-Z]*" ] == ASCIIString["GLY", "LYS"]

  @test ASCIIString[ res.id.name for res in @residues pdb   model "1" chain "A" group "*" residue r"^221[A-Z]?$" ] == ASCIIString["ALA", "LEU"]
  @test ASCIIString[ res.id.name for res in @residues pdbml model "1" chain "A" group "*" residue r"^221[A-Z]?$" ] == ASCIIString["ALA", "LEU"]
end

print("""
1IGY => Insertions have more than one copy of the same amino acid in a single insertion block.
For example, chain B in 1IGY contains a block of four residues inserted at sequence position 82. The block contains Leu-Ser-Ser-Leu.
""")

let code = "1IGY"
  pdb = read(txt(code), PDBFile)
  pdbml = read(xml(code), PDBML)

  @test ASCIIString[ res.id.name for res in @residues pdb   model "1" chain "B" group "*" residue r"^82[A-Z]?$" ] == ["LEU", "SER", "SER", "LEU"]
  @test ASCIIString[ res.id.name for res in @residues pdbml model "1" chain "B" group "*" residue r"^82[A-Z]?$" ] == ["LEU", "SER", "SER", "LEU"]

  @test sum( [ res.id.group for res in @residues pdb model "1" chain "D" group "*" residue "*" ] .== "HETATM" ) == length( @residues pdb model "1" chain "D" group "HETATM" residue "*" )
  @test sum( [ res.id.group for res in @residues pdb model "1" chain "D" group "*" residue "*" ] .== "ATOM" ) == length( @residues pdb model "1" chain "D" group "ATOM" residue "*" )
end

print("""
1HAG => Chain E begins with 1H, 1G, 1F, ... 1A, then 1 (in reverse alphabetic order)
""")

let code = "1HAG"
  pdb = read(txt(code), PDBFile)
  pdbml = read(xml(code), PDBML)

  @test unique(ASCIIString[ res.id.chain for res in pdb   ]) == ["E", "I"]
  @test unique(ASCIIString[ res.id.chain for res in pdbml ]) == ["A", "B", "C", "D", "E"]

  # The chain E of PDB is the chain A of PDBML
  @test ASCIIString[ res.id.number for res in @residues pdb   model "1" chain "E" group "*" residue r"^1[A-Z]?$" ] == ASCIIString[ string(1, code) for code in vcat(collect('H':-1:'A'), "") ]
  @test ASCIIString[ res.id.number for res in @residues pdbml model "1" chain "A" group "*" residue r"^1[A-Z]?$" ] == ASCIIString[ string(1, code) for code in vcat(collect('H':-1:'A'), "") ]
end

print("""
1NSA => Contains a single (unnamed) protein chain with sequence 7A-95A that continues 4-308.
""")

let code = "1NSA"
  pdb = read(txt(code), PDBFile)
  pdbml = read(xml(code), PDBML)

  # Single "A" chain for PDB (auth_asym_id in PDBML)
  @test unique(ASCIIString[ res.id.chain for res in pdb  ]) == ASCIIString[ "A" ]
  # But 'A':'F' chains for PDBML (label_asym_id)
  @test unique(ASCIIString[ res.id.chain for res in pdbml]) == ASCIIString[ string(chain) for chain in 'A':'F' ]

  ind = findobjects(pdbml, Is(:number, "95A"))[1]
  @test pdbml[ind + 1].id.number == "4"

  ind = findobjects(pdb, Is(:number, "95A"))[1]
  @test pdb[ind + 1].id.number == "4"
end

print("""
1IAO => Contains in chain B (in this order) 1S, 2S, 323P-334P, 6-94, 94A, 95-188, 1T, 2T
""")

let code = "1IAO"
  pdb = read(txt(code), PDBFile)
  pdbml = read(xml(code), PDBML)

  pdb_B   = @residues pdb   model "1" chain "B" group "*" residue "*"
  pdbml_B = @residues pdbml model "1" chain "B" group "*" residue "*"

  @test pdb_B[   findobjects(pdb_B,   Is(:number, "2S"))[1]   + 1].id.number == "323P"
  @test pdb_B[   findobjects(pdb_B,   Is(:number, "334P"))[1] + 1].id.number == "6"
  @test pdb_B[   findobjects(pdb_B,   Is(:number, "94"))[1]   + 1].id.number == "94A"
  @test pdb_B[   findobjects(pdb_B,   Is(:number, "94A"))[1]  + 1].id.number == "95"
  @test pdb_B[   findobjects(pdb_B,   Is(:number, "188"))[1]  + 1].id.number == "1T"
  @test pdb_B[   findobjects(pdb_B,   Is(:number, "1T"))[1]   + 1].id.number == "2T"

  @test pdbml_B[ findobjects(pdbml_B, Is(:number, "2S"))[1]   + 1].id.number == "323P"
  @test pdbml_B[ findobjects(pdbml_B, Is(:number, "334P"))[1] + 1].id.number == "6"
  @test pdbml_B[ findobjects(pdbml_B, Is(:number, "94"))[1]   + 1].id.number == "94A"
  @test pdbml_B[ findobjects(pdbml_B, Is(:number, "94A"))[1]  + 1].id.number == "95"
  @test pdbml_B[ findobjects(pdbml_B, Is(:number, "188"))[1]  + 1].id.number == "1T"
  @test pdbml_B[ findobjects(pdbml_B, Is(:number, "1T"))[1]   + 1].id.number == "2T"
end

print("""

RESTful PDB Interface
=====================
""")

@test getpdbdescription("4HHB")["resolution"] == "1.74"
@test_throws KeyError getpdbdescription("104D")["resolution"] # NMR
