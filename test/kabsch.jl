print("""

Test Kabsch algorithm
=====================
""")

let a = [ 1  1 0
         2  1 0
         1+cos(pi/4) 1-sin(pi/4) 0
         1  0 0 ],
    b = [  1 1 0
         1 2 0
         1+cos(pi/4) 1+sin(pi/4)  0 ], # mean(b,1) 1.2357  1.56904  0.0
    sa = a[1:3, :],
    ma = mean(sa, 1) # mean x, y, z for a[1:3,:]

    # Center b, sa, a
    center!(b)
    sa[:,:] = sa .- ma
    a[:,:]  = a  .- ma # center != 0

    @test_approx_eq_eps mean(b,1) [0.0, 0.0, 0.0] 1e-14
    @test_approx_eq_eps mean(sa,1) [0.0, 0.0, 0.0] 1e-14

    F = kabsch(b, sa) # Reference: b
    R = sa * F
    RR = a * F

    @test_approx_eq R b
    @test_approx_eq RR[1:3,:] b

    @test_approx_eq RR[4,:] .- RR[1,:] [1.0, 0.0, 0.0]
    @test_approx_eq sqrt( (RR[1,1] - RR[4,1])^2 + (RR[1,2] - RR[4,2])^2 ) 1.0

    @test_approx_eq_eps rmsd(R, b) 0.0 1e-14
    @test_approx_eq_eps rmsd(RR[1:3,:], b) 0.0 1e-14
end

let P = [51.65 -1.90 50.07
         50.40 -1.23 50.65
         50.68 -0.04 51.54
         50.22 -0.02 52.85],
    Q = [51.30 -2.99 46.54
         51.09 -1.88 47.58
         52.36 -1.20 48.03
         52.71 -1.18 49.38], # P and Q are from BiomolecularStructures.jl' kabsch tests
    Qdistances = Float64[ sqrt(sumabs2(Q[j,:] .- Q[i,:])) for i in 1:4, j in 1:4 ]

    center!(P)
    center!(Q)

    @test_approx_eq_eps mean(P, 1) zeros(3) 1e-14
    @test_approx_eq_eps mean(Q, 1) zeros(3) 1e-14

    rotationmatrix = kabsch(P, Q)
    rotated = Q * rotationmatrix

    @test_approx_eq rmsd(P, rotated) 0.0030426652601371583

    # Internal distances musn't change
    @test_approx_eq Float64[ sqrt(sumabs2(Q[j,:] .- Q[i,:])) for i in 1:4, j in 1:4 ] Qdistances # Translated
    @test_approx_eq Float64[ sqrt(sumabs2(rotated[j,:] .- rotated[i,:])) for i in 1:4, j in 1:4 ] Qdistances # Translated  and rotated
end

let P = [-1. 0. 0.
         0. 2. 0.
         0. 1. 0.
         0. 1. 1],
    Q = [0. -1. -1.
         0. -1. 0
         0. 0. 0.
         -1. 0. 0.]

    center!(P)
    center!(Q)
    rotated = Q * kabsch(P, Q)
    @test_approx_eq_eps rmsd(P, rotated) 0.695 0.001
end

print("""

Superimpose PDBs and RMSD
-------------------------
""")

let hemoglobin = read(joinpath(pwd(), "data", "2hhb.pdb.gz"), PDBFile, group="ATOM", model="1")

    α1 = @residues hemoglobin model "1" chain "A" group "ATOM" residue "*"
    α2 = @residues hemoglobin model "1" chain "C" group "ATOM" residue "*"
    β1 = @residues hemoglobin model "1" chain "B" group "ATOM" residue "*"
    β2 = @residues hemoglobin model "1" chain "D" group "ATOM" residue "*"

    a1, a2, rα = superimpose(α1, α2)

    @test_approx_eq_eps rα 0.230 0.001 # Bio3D RMSD

    b1, b2, rβ = superimpose(β1, β2)

    @test_approx_eq_eps rβ 0.251 0.001 # Bio3D & Chimera's MatchMaker RMSD

    @test length(a1) == length(α1)
    @test length(a2) == length(α2)
    @test length(b1) == length(β1)
    @test length(b2) == length(β2)

    @test a1 != α1
    @test a2 != α2

    @test b1 != β1
    @test b2 != β2

    ca_a1 = CAmatrix(a1)
    ca_a2 = CAmatrix(a2)
    ca_b1 = CAmatrix(b1)
    ca_b2 = CAmatrix(b2)

    @test PDB._iscentered(ca_a1)
    @test PDB._iscentered(ca_a2)
    @test PDB._iscentered(ca_b1)
    @test PDB._iscentered(ca_b2)

    @test_approx_eq_eps mean(ca_a1,1) zeros(3) 1e-14
    @test_approx_eq_eps mean(ca_a2,1) zeros(3) 1e-14
    @test_approx_eq_eps mean(ca_b1,1) zeros(3) 1e-14
    @test_approx_eq_eps mean(ca_b2,1) zeros(3) 1e-14

    @test rα == rmsd(ca_a1, ca_a2)
    @test rβ == rmsd(ca_b1, ca_b2)
end

