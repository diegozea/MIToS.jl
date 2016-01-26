# using Base.Test
# using MIToS.MSA

print("""

Tests for Raw Multiple Sequence Alignment a.k.a Matrix{Residue}
================================================================
""")

print("""

Parse Raw
---------
""")

let raw = read(joinpath(pwd(), "data", "gaps.txt"), Raw),
  raw_string = """THAYQAIHQV
       THAYQAIHQ-
       THAYQAIH--
       THAYQAI---
       THAYQA----
       THAYQ-----
       THAY------
       THA-------
       TH--------
       T---------
       """

  @test asciisequence(raw, 1)  == "THAYQAIHQV"
  @test asciisequence(raw, 10) == "T---------"

  @test parse(raw_string, Raw) == raw

  print("""

Print Raw
---------
""")

  let io = IOBuffer()
    print(io, raw, Raw)
    @test takebuf_string(io) == raw_string
  end

  print("""

Test %
------
""")

  @test gappercentage(raw) == 0.45
  @test gappercentage(raw, 1) == [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
  @test gappercentage(raw, 2) == [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]

  @test residuepercentage(raw) == 0.55
  @test residuepercentage(raw, 1) == [1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1]
  @test residuepercentage(raw, 2) == [1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1]

  @test coverage(raw) == [1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1]

  @test columngappercentage(raw) == [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]

  print("""

Reference and Gapstrip
----------------------
""")

  let gs = gapstrip(raw,  coveragelimit=0.5, gaplimit=0.5)
    @test getsequence(gs, 1) == res"THAYQAIH"
    @test ncolumns(gs) == 8
    @test nsequences(gs) == 6
  end

  let ref = setreference!(copy(raw), 2)
    @test getsequence(ref, 1) == res"THAYQAIHQ-"
    ref = adjustreference(ref)
    @test getsequence(ref, 1) == res"THAYQAIHQ"
  end

end
