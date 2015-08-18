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

let map = siftsmapping(sifts_file, 'A', "2vqc", db="PDB")
  @test_throws KeyError map[9]  # Missing
  @test_throws KeyError map[80] # Missing
  @test_throws KeyError map[1]  # Missing
  @test map[10] == 4
  @test map[79] == 73
end

let map = siftsmapping(sifts_file, 'A', "2vqc", db="PDB", check_observed=false)
  @test map[9] == 3   # Missing
  @test map[80] == 74 # Missing
  @test map[1] == -5  # Missing
  @test map[10] == 4
  @test map[79] == 73
end
