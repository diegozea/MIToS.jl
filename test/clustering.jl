# using Base.Test
# using MIToS.Clustering
# using MIToS.MSA

println("""

Clustering
==========
""")

println("""

percentidentity
---------------
""")

print("""
-> Float64
""")

@test_throws ErrorException percentidentity(res"AH", res"AGH")

@test percentidentity(res"AH", res"AH") == 1.
@test percentidentity(res"AH", res"AG") == .5
@test percentidentity(res"AH", res"RG") == 0.
@test percentidentity(res"AH-", res"AG-") == .5
@test percentidentity(res"A--", res"AG-") == .5

print("""
-> Bool
""")

@test percentidentity(res"A--", res"AG-", 0.4)
@test !percentidentity(res"A--", res"AG-", 0.6)

print("""

Hobohm I
--------
""")

# DAWAEE
# DAWAEF  83.3
# DAWAED  83.3
# DAYCMD  33.3
# DAYCMT  33.3  83.3
# DAYCMT  33.3  83.3

let fasta = read(joinpath(pwd(), "data", "Gaoetal2011.fasta"), FASTA)
  clusters = hobohmI(fasta, 0.62)

  @test getnclusters(clusters) == 2
  @test nsequences(clusters) == 6
  @test getweight(clusters, 1) == 1/3
  @test getweight(clusters, 6) == 1/3
end

print("""

percentidentity on a MSA
------------------------
""")

let fasta = read(joinpath(pwd(), "data", "Gaoetal2011.fasta"), FASTA)
  id = percentidentity(fasta)

  @test id[1,1] == 100.0
  @test_approx_eq_eps id[1,2] 83.33 0.01
  @test_approx_eq_eps id[1,3] 83.33 0.01
  @test_approx_eq_eps id[3,1] 83.33 0.01
  @test_approx_eq_eps id[1,6] 33.33 0.01
  @test_approx_eq_eps id[4,5] 83.33 0.01
  @test id[5,6] == 100.0

  @test maximum(id) == 100.0
  @test_approx_eq_eps minimum(id) 33.33 0.01

  print("""
  SequenceIdentityMatrix
  """)
  @test isapprox(sum(id[:,3]), 100.0 + 50.0 + 2*(200/6) + 2*(500/6))
  id[4,3] = 80
  @test isapprox(sum(id[:,3]), 100.0 + 80.0 + 2*(200/6) + 2*(500/6))
end

let aln = read(joinpath(pwd(), "data", "gaps.txt"), Raw)
  id = percentidentity(aln)

  @test id[1,1] == 100.0
  @test id[1,2] == 90.0
  @test id[1,3] == 80.0
  @test_approx_eq id[2,3] 800/9
end
