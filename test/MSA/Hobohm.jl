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
    @test nsequences(clusters) == 6
    @test getweight(clusters, 1) == 1/3
    @test getweight(clusters, 6) == 1/3
end
