# using Base.Test
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

@test percentidentity(res"AH", res"AH") == 100.
@test percentidentity(res"AH", res"AG") == 50.
@test percentidentity(res"AH", res"RG") == 0.
@test percentidentity(res"AH-", res"AG-") == 50.
@test percentidentity(res"A--", res"AG-") == 50.

print("""
-> Bool
""")

@test percentidentity(res"A--", res"AG-", 40.)
@test !percentidentity(res"A--", res"AG-", 60)

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
  clusters = hobohmI(fasta, 62)

  @test nclusters(clusters) == 2
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

print("""

Test meanpercentidentity
------------------------
""")

let msa = vcat( transpose(res"--GGG-"),
                transpose(res"---GGG") )
#                  identities 000110 sum 2
#            aligned residues 001111 sum 4

    @test percentidentity(msa)[1, 2] == 50.0 # 2 / 4
    @test meanpercentidentity(msa)   == 50.0
end

let msa    = rand(Residue, 400, 2),
    msa300 = msa[1:300, :]

    @test mean(percentidentity(msa300).list) == meanpercentidentity(msa300)
    @test_approx_eq_eps mean(percentidentity(msa).list) meanpercentidentity(msa) 0.5
    @test mean(percentidentity(msa).list) == meanpercentidentity(msa, exact=true)
end

print("""

Test percentsimilarity
======================
""")

@test percentsimilarity(res"AH", res"IM") == 50.0

print("""
Test with gaps
""")

@test percentsimilarity(res"-AH", res"--H") == 50.0
@test percentsimilarity(res"-AH", res"-AH") == 100.0
@test percentsimilarity(res"-AH", res"-AH") == percentsimilarity(res"AH", res"AH")
@test percentsimilarity(res"-AH", res"-AH") == percentsimilarity(res"--AH", res"--AH")

print("""
Using SMS's "Ident and Sim" residue groups
""")

let fasta = read(joinpath(pwd(), "data", "Gaoetal2011.fasta"), FASTA),
    sim = percentsimilarity(fasta, [1, 4, 5, 5, 7, 5, 5, 1, 4, 1, 1, 4, 7, 2, 6, 3, 3, 2, 2, 1])

  @test eltype(sim) == Float64

  @test sim[1,1] == 100.0

  @test_approx_eq_eps sim[1,2]  83.33 0.01
  @test_approx_eq_eps sim[1,3] 100.00 0.01
  @test_approx_eq_eps sim[1,4]  66.67 0.01
  @test_approx_eq_eps sim[1,5]  50.00 0.01
  @test_approx_eq_eps sim[1,6]  50.00 0.01
  @test_approx_eq_eps sim[2,3]  83.33 0.01
  @test_approx_eq_eps sim[2,4]  50.00 0.01
  @test_approx_eq_eps sim[2,5]  50.00 0.01
  @test_approx_eq_eps sim[2,6]  50.00 0.01
  @test_approx_eq_eps sim[3,4]  66.67 0.01
  @test_approx_eq_eps sim[3,5]  50.00 0.01
  @test_approx_eq_eps sim[3,6]  50.00 0.01
  @test_approx_eq_eps sim[4,5]  83.33 0.01
  @test_approx_eq_eps sim[4,6]  83.33 0.01
  @test_approx_eq_eps sim[5,6] 100.00 0.01

end

print("""
Using Bio3D's (2.2) seqidentity residue groups
""")

let fasta = read(joinpath(pwd(), "data", "Gaoetal2011.fasta"), FASTA),
    sim = percentsimilarity(fasta, [1, 5, 7, 6, 9, 7, 6, 1, 5, 2, 2, 5, 2, 3, 8, 4, 4, 3, 3, 2], out=Float16),
    bio3d = [ 1.000 0.833 1.000 0.667 0.500 0.500
              0.833 1.000 0.833 0.500 0.500 0.500
              1.000 0.833 1.000 0.667 0.500 0.500
              0.667 0.500 0.667 1.000 0.833 0.833
              0.500 0.500 0.500 0.833 1.000 1.000
              0.500 0.500 0.500 0.833 1.000 1.000 ] .* 100.00

    @test eltype(sim) == Float16

    for i in 1:6
        for j in 1:6
            @test_approx_eq_eps sim[i,j] bio3d[i,j] 0.1
        end
    end

end
