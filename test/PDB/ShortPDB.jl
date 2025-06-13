@testset "Read PDB with short atom lines" begin
    pdbfile = joinpath(DATA, "short.pdb")
    pdb = read_file(pdbfile, PDBFile)
    @test length(pdb) == 1
    atoms = pdb[1].atoms
    @test length(atoms) == 8
    @test atoms[1].occupancy == 1.0
    @test atoms[1].B == "0.0"
    @test atoms[2].occupancy == 1.0
    @test atoms[2].B == "0.0"
    @test atoms[3].occupancy == 1.0
    @test atoms[3].B == "13.79"
    @test atoms[4].occupancy == 1.0
    @test atoms[4].B == "0.0"
    @test atoms[5].occupancy == 0.50
    @test atoms[5].B == "0.0"
    @test atoms[6].occupancy == 1.0
    @test atoms[6].B == "12.34"
    @test atoms[7].occupancy == 0.70
    @test atoms[7].B == "0.0"
    @test atoms[8].occupancy == 1.0
    @test atoms[8].B == "0.0"
end
