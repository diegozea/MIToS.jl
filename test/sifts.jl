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

# 1SSX => 15A 15B
# 1AS5 => NMR
# 1CBN => Rotamers: Thr1 (+Calpha), Thr2, Ile7, Val8, Arg10
# 1DPO => When inserted residues lack insertion letters (probably illegal according to the PDB format specification), PE fails to show the inserted residues in Sequences/Seq3D. For example, the single unnamed chain in 1DPO contains insertions at postions 184 (Gly, Phe), 188 (Gly, Lys), and 221 (Ala, Leu) but no insertion letters. Clicking on Gly184 in Seq3D highlights both Gly184 and Phe184. 
# 1IGY => When insertions have more than one copy of the same amino acid (or nucleotide) in a single insertion block, clicking one copy in Seq3D highlights all copies within the same insertion block. For example, chain B in 1IGY contains a block of four residues inserted at sequence position 82. The block contains Leu-Ser-Ser-Leu. Clicking on either Leu highlights both of them in the molecular view, and similarly for either Ser. This results from the inability of Chime to use insertion codes in selecting atoms. For example, Chime cannot distinguish Leu82 from Leu82C. 
# For example, chain E in 1HAG begins with 1H, 1G, 1F, ... 1A, then 1 (in reverse alphabetic order)
# Sequence numbers in decreasing order within a single chain cause PE's Sequences report to contain errors, and cause Seq3D to fail to work properly. This is due to design flaws in PE. Because these cases are rare, the effort to remedy these flaws doesn't seem worthwhile. An example is 1NSA, which contains a single (unnamed) protein chain with sequence 7A-95A that continues 4-308. Another example is 1IAO, which contains in chain B (in this order) 1S, 323P-334P, 6-94, 94A, 95-188, 1T, 2T. PE detects such situations and alerts the user that the Sequences/Seq3D displays are garbled.
# When a single PDB file contains both an unnamed and named chains, PE's mechanisms that depend upon chain names fail. This is probably quite rare, but occurs in 4CPA, where QuickViews's SELECT menu lists the unnamed chain, but clicking on it selects all chains. Also, in Seq3D, clicking on D16 in the unnamed chain highlights D16 there and in chain I (which also happens to have D at position 16).ZZ
