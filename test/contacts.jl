using Base.Test
using MIToS.PDB

txt(code) = joinpath(pwd(), "data", string(uppercase(code), ".pdb"))
xml(code) = joinpath(pwd(), "data", string(uppercase(code), ".xml"))

print("""

Tests for Contacts on PDB
=========================
""")

print("""
1IGY => Test it using data from http://www-cryst.bioc.cam.ac.uk/~richard/piccolo/piccolo.php?PDB=1IGY (28/Sep/2015)
""")

let code = "1IGY"
  pdb = read(txt(code), PDBFile)
  # pdbml = read(xml(code), PDBML)

  # Contacts: Chain A (4 residues) y Chain D (3 residues)
  C1 = @residuesdict pdb model "1" chain "A" group "ATOM" residue "*"
  C2 = @residuesdict pdb model "1" chain "D" group "ATOM" residue "*"

  TRUE = [ ("126",	"311"), ("183",	"309"), ("183",	"310"),
           ("184",	"309"), ("187",	"309"), ("187",	"310"),
           ("187",	"311"), ("187",	"312"), ("187",	"319"),
           ("210",	"237"), ("211",	"312"), ("212",	"237"),
           ("213",	"237"), ("213",	"312"), ("213",	"313"),
           ("214",	"237") ]

  for (resnum1, resnum2) in TRUE
    res1 = C1[resnum1]
    res2 = C2[resnum2]

    @test contact(res1, res2, 6.5)

    if resnum2 == "309" && (resnum1 == "184" || resnum1 == "187")
      @test ionic(res1, res2)
    else
      @test !ionic(res1, res2)
    end

    if (resnum1 == "211" && resnum2 == "312") || (resnum1 == "212" && resnum2 == "237")
      @test vanderwaals(res1, res2)
    else
      @test !vanderwaals(res1, res2)
    end

    if resnum1 == "212" && resnum2 == "237"
      @test hydrophobic(res1, res2)
    else
      @test !hydrophobic(res1, res2)
    end
    
    if resnum1 == "211" && resnum2 == "312"
      @test vanderwaalsclash(res1, res2)
    else
      @test !vanderwaalsclash(res1, res2)
    end

    @test !aromaticsulphur(res1, res2)
    @test !pication(res1, res2)
    @test !disulphide(res1, res2)
    @test !aromatic(res1, res2)
    @test !hydrogenbond(res1, res2)
    @test !hydrogenbond(res1, res2)
    @test !covalent(res1, res2)
    
  end

end
