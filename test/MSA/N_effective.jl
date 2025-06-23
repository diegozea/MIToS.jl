@testset "n_effective" begin
    msa = reshape([res"A"; res"C"; res"D"], 3, 1)
    @test n_effective(msa, 0.0) == 1.0
    @test n_effective(msa, 100.0) == 3.0
end
