@testset "Kabsch algorithm" begin

    @testset "Test I" begin
        a = [1  1 0
             2  1 0
             1+cos(pi/4) 1-sin(pi/4) 0
             1  0 0 ]
        b = [1 1 0
             1 2 0
             1+cos(pi/4) 1+sin(pi/4)  0 ] # mean(b,1) 1.2357  1.56904  0.0
        sa = a[1:3, :]
        ma = mean(sa, dims=1) # mean x, y, z for a[1:3,:]

        # Center b, sa, a
        center!(b)
        sa[:,:] = sa .- ma
        a[:,:]  = a  .- ma # center != 0

        @test isapprox(vec(mean(b, dims=1)), [0.0, 0.0, 0.0], atol=1e-13)
        @test isapprox(vec(mean(sa, dims=1)), [0.0, 0.0, 0.0], atol=1e-13)

        F = kabsch(b, sa) # Reference: b
        R = sa * F
        RR = a * F

        @test R ≈ b
        @test RR[1:3,:] ≈ b

        @test RR[4,:] .- RR[1,:] ≈ [1.0, 0.0, 0.0]
        @test sqrt((RR[1,1] - RR[4,1])^2 + (RR[1,2] - RR[4,2])^2) ≈ 1.0

        @test isapprox(PDB.rmsd(R, b), 0.0, atol=1e-13)
        @test isapprox(PDB.rmsd(RR[1:3,:], b), 0.0, atol=1e-13)
    end

    @testset "BiomolecularStructures' test" begin

        P = [51.65 -1.90 50.07
             50.40 -1.23 50.65
             50.68 -0.04 51.54
             50.22 -0.02 52.85]
        Q = [51.30 -2.99 46.54
             51.09 -1.88 47.58
             52.36 -1.20 48.03
             52.71 -1.18 49.38] # P and Q are from BiomolecularStructures.jl' kabsch tests
        Qdistances = Float64[ sqrt(sum(abs2, Q[j,:] .- Q[i,:])) for i in 1:4, j in 1:4 ]

        center!(P)
        center!(Q)

        @test isapprox(vec(mean(P, dims=1)), zeros(3), atol=1e-13)
        @test isapprox(vec(mean(Q, dims=1)), zeros(3), atol=1e-13)

        rotationmatrix = kabsch(P, Q)
        rotated = Q * rotationmatrix

        @test PDB.rmsd(P, rotated) ≈ 0.0030426652601371583

        # Internal distances mustn't change
        @test Float64[sqrt(sum(abs2, Q[j,:] .- Q[i,:])) for i in 1:4, j in 1:4] ≈ Qdistances # Translated
        @test Float64[sqrt(sum(abs2, rotated[j,:] .- rotated[i,:])) for i in 1:4, j in 1:4] ≈ Qdistances # Translated  and rotated
    end

    @testset "Test II" begin
        P = [-1. 0. 0.
             0. 2. 0.
             0. 1. 0.
             0. 1. 1]
        Q = [0. -1. -1.
             0. -1. 0
             0. 0. 0.
             -1. 0. 0.]

        center!(P)
        center!(Q)
        rotated = Q * kabsch(P, Q)
        @test isapprox(PDB.rmsd(P, rotated), 0.695, atol=0.001)
    end

    @testset "RMSF" begin
        A = [ -1. -1. -1.
              -1.  0.  1. ]
        B = [  1.  1.  1.
              -1.  0.  1. ]
        w = [ .25 , .75 ]

        @test_throws ArgumentError mean_coordinates(Matrix{Float64}[A[1:1,:], B])
        @test_throws ArgumentError mean_coordinates(Matrix{Float64}[A[:,1:2], B])
        @test_throws ArgumentError mean_coordinates(Matrix{Float64}[A])
        @test_throws ArgumentError mean_coordinates(Matrix{Float64}[A,  B'])
        @test_throws ArgumentError mean_coordinates(Matrix{Float64}[A', B'])

        @test mean_coordinates(Matrix{Float64}[A, B]) == [ 0. 0. 0.
                                                          -1. 0. 1. ]

        @test mean_coordinates(Matrix{Float64}[A, B], w) == [ 0.5 0.5 0.5
                                                             -1.  0.  1. ]

        @test mean_coordinates(Matrix{Float64}[A, B], reverse(w)) == [ -0.5 -0.5 -0.5
                                                                       -1.   0.   1. ]

        @test rmsf(Matrix{Float64}[A, B]) == [sqrt(3), 0.]
        @test rmsf(Matrix{Float64}[A, B], w) == [sqrt(0.25*6.75 + 0.75*0.75), 0.]
    end

    @testset "Superimpose PDBs" begin

        hemoglobin = read(joinpath(pwd(),"data","2hhb.pdb.gz"),PDBFile,group="ATOM",model="1")

        α1 = @residues hemoglobin model "1" chain "A" group "ATOM" residue All
        α2 = @residues hemoglobin model "1" chain "C" group "ATOM" residue All
        β1 = @residues hemoglobin model "1" chain "B" group "ATOM" residue All
        β2 = @residues hemoglobin model "1" chain "D" group "ATOM" residue All

        a1, a2, rα = superimpose(α1, α2)

        @test isapprox(rα, 0.230, atol=0.001) # Bio3D RMSD

        b1, b2, rβ = superimpose(β1, β2)

        @test isapprox(rβ, 0.251, atol=0.001) # Bio3D & Chimera's MatchMaker RMSD

        @test length(a1) == length(α1)
        @test length(a2) == length(α2)
        @test length(b1) == length(β1)
        @test length(b2) == length(β2)

        @test a1 != α1
        @test a2 != α2

        @test b1 != β1
        @test b2 != β2

        @testset "Vector{PDBResidue}" begin

            @test isapprox(PDB.rmsd(α1, α2), 0.230, atol=0.001)
            @test isapprox(PDB.rmsd(β1, β2), 0.251, atol=0.001)

            @test isapprox(PDB.rmsd(a1, a2, superimposed=true), 0.230, atol=0.001)
            @test isapprox(PDB.rmsd(b1, b2, superimposed=true), 0.251, atol=0.001)

            ca_a1 = CAmatrix(a1)
            ca_a2 = CAmatrix(a2)
            ca_b1 = CAmatrix(b1)
            ca_b2 = CAmatrix(b2)

            @testset "Centered" begin
                @test isapprox(mean(ca_a1, dims=1), zeros(1,3), atol=1e-13)
                @test isapprox(mean(ca_a2, dims=1), zeros(1,3), atol=1e-13)
                @test isapprox(mean(ca_b1, dims=1), zeros(1,3), atol=1e-13)
                @test isapprox(mean(ca_b2, dims=1), zeros(1,3), atol=1e-13)

                @test PDB._iscentered(ca_a1)
                @test PDB._iscentered(ca_a2)
                @test PDB._iscentered(ca_b1)
                @test PDB._iscentered(ca_b2)
            end

            @testset "RMSD" begin

                @test rα == PDB.rmsd(ca_a1, ca_a2)
                @test rβ == PDB.rmsd(ca_b1, ca_b2)
            end
        end

        @testset "coordinatesmatrix and centered..." begin

            @test coordinatesmatrix(centeredresidues(α1)) ≈ coordinatesmatrix(a1)
            @test coordinatesmatrix(centeredresidues(β1)) ≈ coordinatesmatrix(b1)

            @test centeredcoordinates(α1) ≈ coordinatesmatrix(a1)
            @test centeredcoordinates(β1) ≈ coordinatesmatrix(b1)
        end

        @testset "Superimpose with centered PDBs" begin

            b12, b22, rβ2 = superimpose(b1, b2)
            @test rβ2 ≈ rβ
        end

        @testset "RMSF" begin

            @test rmsf(Vector{PDBResidue}[α1,α2]) == rmsf(Vector{PDBResidue}[α1,α2],[.5,.5])
            @test rmsf(Vector{PDBResidue}[β1,β2]) == rmsf(Vector{PDBResidue}[β1,β2],[.5,.5])
        end
    end

    @testset "PDBResidue without alpha-carbon" begin
        small_2WEL = read(joinpath(pwd(), "data", "2WEL_D_region.pdb"), PDBFile)
        small_6BAB = read(joinpath(pwd(), "data", "6BAB_D_region.pdb"), PDBFile)

        aln_2WEL, aln_6BAB, RMSD = superimpose(small_2WEL, small_6BAB)
        @test length(aln_2WEL) == 3
        @test length(aln_6BAB) == 3
        @test small_2WEL[2].atoms[1].coordinates ≉ aln_2WEL[2].atoms[1].coordinates
        @test small_6BAB[2].atoms[1].coordinates ≉ aln_6BAB[2].atoms[1].coordinates
        @test RMSD < 1.0e-14 # e.g. 6.9e-15, 3.7e-15 & 1.6e-15

        @testset "PDB without CA" begin

            filter!(atom -> atom.atom != "CA", small_2WEL[1].atoms)
            filter!(atom -> atom.atom != "CA", small_6BAB[3].atoms)
            @test_throws ArgumentError superimpose(small_2WEL, small_6BAB)
        end
    end
end
