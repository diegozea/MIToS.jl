@testset "Clusters" begin

    for int = 1:100
        @test getweight(NoClustering(), int) == 1.0
    end
end

@testset "Hobohm I" begin
    cluster_weights(cl) = [getweight(cl, i) for i = 1:nelements(cl)]

    # DAWAEE
    # DAWAEF  83.3
    # DAWAED  83.3
    # DAYCMD  33.3
    # DAYCMT  33.3  83.3
    # DAYCMT  33.3  83.3

    fasta = read_file(joinpath(DATA, "Gaoetal2011.fasta"), FASTA)
    clusters = hobohmI(fasta, 62)

    @test nclusters(clusters) == 2
    @test nelements(clusters) == 6
    @test getweight(clusters, 1) == 1 / 3
    @test getweight(clusters, 6) == 1 / 3

    @testset "Serial and threaded are exact" begin
        for (file, format, threshold) in
            [("Gaoetal2011.fasta", FASTA, 62), ("gaps.txt", Raw, 62)]
            msa = read_file(joinpath(DATA, file), format)
            serial = hobohmI(msa, threshold; threads = false)
            threaded = hobohmI(msa, threshold; threads = true)

            @test serial == threaded
            @test assignments(serial) == assignments(threaded)
            @test counts(serial) == counts(threaded)
            @test cluster_weights(serial) == cluster_weights(threaded)
        end
    end

    @testset "Clusters getters" begin

        @test getweight(clusters) == clusters.weights
        @test assignments(clusters) == clusters.clusters
        @test counts(clusters) == clusters.clustersize
    end

    @testset "Convert to Clusters" begin

        @test convert(Clusters, clusters) == clusters

        distance = convert(Matrix{Float64}, 100.0 .- percentidentity(fasta))
        cr = Clustering.dbscan(distance, 38.0, metric = nothing, min_neighbors = 2)
        @test convert(Clusters, cr) == clusters
    end

    @testset "Do-block predicate" begin
        clusters_do = hobohmI(fasta, 62) do s1, s2, thr
            percentidentity(s1, s2, thr)
        end
        @test clusters_do == clusters

        clusters_do_serial = hobohmI(fasta, 62; threads = false) do s1, s2, thr
            percentidentity(s1, s2, thr)
        end
        clusters_do_threaded = hobohmI(fasta, 62; threads = true) do s1, s2, thr
            percentidentity(s1, s2, thr)
        end
        @test clusters_do == clusters_do_serial
        @test clusters_do_serial == clusters_do_threaded
    end

    @testset "Explicit percentidentity keeps threaded default" begin
        clusters_percentidentity = hobohmI(percentidentity, fasta, 62)
        clusters_percentidentity_serial =
            hobohmI(percentidentity, fasta, 62; threads = false)
        clusters_percentidentity_threaded =
            hobohmI(percentidentity, fasta, 62; threads = true)

        @test clusters_percentidentity == clusters
        @test clusters_percentidentity == clusters_percentidentity_threaded
        @test clusters_percentidentity_serial == clusters_percentidentity_threaded
    end

    @testset "Vector input" begin
        seqs = getresiduesequences(fasta)
        clusters_vec = hobohmI(percentidentity, seqs, 62)
        @test clusters_vec == clusters

        clusters_vec_serial = hobohmI(percentidentity, seqs, 62; threads = false)
        clusters_vec_threaded = hobohmI(percentidentity, seqs, 62; threads = true)
        @test clusters_vec == clusters_vec_threaded
        @test clusters_vec_serial == clusters_vec_threaded
    end

    @testset "Custom predicates stay serial by default" begin
        seqs = [res"AAAA" for _ = 1:128]

        function serial_only_percentidentity(s1, s2, thr)
            Threads.threadid() == 1 || error("predicate ran on a worker thread")
            percentidentity(s1, s2, thr)
        end

        clusters_default = hobohmI(serial_only_percentidentity, seqs, 100)
        clusters_serial = hobohmI(serial_only_percentidentity, seqs, 100; threads = false)

        @test clusters_default == clusters_serial
        @test counts(clusters_default) == [128]
    end
end
