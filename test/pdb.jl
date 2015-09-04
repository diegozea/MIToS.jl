using Base.Test
using MIToS.PDB

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

print("""
1AS5 => NMR
""")

print("""
1DPO => Inserted residues lack insertion letters
Single unnamed chain in 1DPO contains insertions at postions 184 (Gly, Phe), 188 (Gly, Lys), and 221 (Ala, Leu) but no insertion letters.
""")

print("""
1IGY => Insertions have more than one copy of the same amino acid in a single insertion block.
For example, chain B in 1IGY contains a block of four residues inserted at sequence position 82. The block contains Leu-Ser-Ser-Leu.
""")

print("""
1HAG => Chain E begins with 1H, 1G, 1F, ... 1A, then 1 (in reverse alphabetic order)
""")

print("""
1NSA => Contains a single (unnamed) protein chain with sequence 7A-95A that continues 4-308.
""")

print("""
1IAO => Contains in chain B (in this order) 1S, 323P-334P, 6-94, 94A, 95-188, 1T, 2T
""")
