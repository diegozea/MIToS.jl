@testset "Three letters name" begin
    
    @test three2residue("ALA") == Residue('A')
    @test three2residue("Ala") == Residue('A')
    @test three2residue("ala") == Residue('A')

    @test_throws ErrorException three2residue("AL")
    @test_throws ErrorException three2residue("Alanine")

    @test residue2three(Residue('A')) == "ALA"

    @test_throws ErrorException residue2three(GAP)
    @test_throws ErrorException residue2three(INV)

    @test three2residue("XAA") == XAA
    @test three2residue("UNK") == XAA
end
