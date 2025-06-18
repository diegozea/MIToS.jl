@testset "MSA Plots recipe" begin
    using RecipesBase
    msa = read_file(joinpath(DATA, "simple.fasta"), FASTA)
    result = RecipesBase.apply_recipe(Dict{Symbol,Any}(), msa)
    @test length(result) == 1
    data = result[1]
    @test data.plotattributes[:seriestype] == :heatmap
    @test data.plotattributes[:yflip]
    @test data.plotattributes[:grid] == false
    @test data.args[1] == 1:ncolumns(msa)
    @test data.args[2] == sequencenames(msa)
    @test data.args[3] == map(string, getresidues(msa))
end

