using Base.Test
using MIToS.Clustering
using MIToS.MSA

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
# DAWAEF  100
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
