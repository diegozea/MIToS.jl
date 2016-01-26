# using Base.Test
# using MIToS.MSA

let alphabet = "ARNDCQEGHILKMFPSTWYV-",
  msa =  ["DAWAEF",
          "DAWAED",
          "DAYCMD" ],
  badmsa=["DAWAEF",
          "DAWAED",
          "DAM"    ]

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
  @test string(Residue(alphabet)) == alphabet # "AR..." instead of the standar "[A,R,..."

  @test convert(Matrix{Residue}, msa) == Residue['D' 'A' 'W' 'A' 'E' 'F'
                                                 'D' 'A' 'W' 'A' 'E' 'D'
                                                 'D' 'A' 'Y' 'C' 'M' 'D']

  @test_throws ErrorException convert(Matrix{Residue}, badmsa)

  print("""
Ambiguous or not standard residues (are gaps on MIToS)
""")

  for res in Char['X', # Any amino acid
                  'B', # Aspartic acid or Asparagine
                  'Z', # Glutamine or Glutamic acid
                  'O', # Pyrrolysine
                  'U', # Selenocysteine
                  'J'] # Leucine or Isoleucine

    @test Residue(res) == GAP

  end

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

  print("""
Random
""")
  for i in 1:8000
    @test rand(Residue) != GAP
    @test typeof( rand(Residue) ) == Residue
  end
  @test size(rand(Residue, 20)) == (20,)
  @test typeof(rand(Residue, 20)) == Array{Residue, 1}
  @test size(rand(Residue, 20,30)) == (20,30)
  @test typeof(rand(Residue, 20, 30)) == Array{Residue, 2}

end
