@testset "Plot recipe" begin
    using MIToS.PDB
    using RecipesBase

    residue = PDBResidue(
        PDBResidueIdentifier("1", "1", "ALA", "ATOM", "1", "A"),
        [PDBAtom(Coordinates(1.0, 2.0, 3.0), "CA", "C", 1.0, "0", "", "")],
    )
    recipe = RecipesBase.apply_recipe(Dict{Symbol,Any}(), [residue])
    @test length(recipe) == 1
    rd = recipe[1]
    @test rd.plotattributes[:group] == ["A"]
    @test rd.args == ([1.0], [2.0], [3.0])
end
