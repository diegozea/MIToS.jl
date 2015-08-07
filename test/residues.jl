using Base.Test
using MIToS.MSA

const alphabet = "ARNDCQEGHILKMFPSTWYV-"

print("""

Test Residue
============
""")

print("""
Convert
""")
@test Int[ Residue(char) for char in alphabet] == Int[ i for i in 1:21 ]
@test Int[ Residue(char) for char in lowercase("ARNDCQEGHILKMFPSTWYV.")] == Int[ 21 for i in 1:21 ]

print("""
GAP
""")
@test Int(GAP) == 21
@test Int8(GAP) == Int8(21)
@test Residue('-') == GAP
@test Char(GAP) == '-'
@test UInt8(GAP) == UInt8('-')

print("""
Convert to and from ASCIIString
""")
@test res"ARNDCQEGHILKMFPSTWYV-" == Residue[ Residue(char) for char in alphabet]
@test Residue(alphabet) == Residue[ Residue(char) for char in alphabet]
@test ascii(Residue(alphabet)) == alphabet

const msa =  ["DAWAEF",
              "DAWAED",
              "DAYCMD" ]

const badmsa=["DAWAEF",
              "DAWAED",
              "DAM"    ]

@test convert(Matrix{Residue}, msa) == Residue['D' 'A' 'W' 'A' 'E' 'F'
                                              'D' 'A' 'W' 'A' 'E' 'D'
                                              'D' 'A' 'Y' 'C' 'M' 'D']

@test_throws ErrorException convert(Matrix{Residue}, badmsa)

print("""
Comparisons
""")

for res in Residue(alphabet)
  @test res == res
end

for i in 1:21, j in 1:21
  if j != i
    @test Residue(i) != Residue(j)
  end
end

