@testset "Clusters" begin

    for int in 1:100
        @test getweight(NoClustering(), int) == 1.0
    end
end

@testset "Hobohm I" begin

    # DAWAEE
    # DAWAEF  83.3
    # DAWAED  83.3
    # DAYCMD  33.3
    # DAYCMT  33.3  83.3
    # DAYCMT  33.3  83.3

    fasta = read(joinpath(pwd(), "data", "Gaoetal2011.fasta"), FASTA)
    clusters = hobohmI(fasta, 62)

    @test nclusters(clusters) == 2
    @test nelements(clusters) == 6
    @test getweight(clusters, 1) == 1/3
    @test getweight(clusters, 6) == 1/3

    @testset "Clusters getters" begin

        @test getweight(clusters) == clusters.weights
        @test assignments(clusters) == clusters.clusters
        @test counts(clusters) == clusters.clustersize
    end

    @testset "Convert to Clusters" begin

        @test convert(Clusters, clusters) == clusters

        distance = convert(Matrix{Float64}, 100.0 .- percentidentity(fasta))
        cr = Clustering.dbscan(distance, 38.0, 2)
        @test convert(Clusters, cr) == clusters
    end
end
