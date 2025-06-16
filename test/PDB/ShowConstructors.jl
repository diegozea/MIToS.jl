@testset "Show constructors" begin
    pdbfile = joinpath(DATA, "short.pdb")
    res = read_file(pdbfile, PDBFile)[1]

    id_str = sprint(show, res.id)
    id = eval(Meta.parse(id_str))
    @test id == res.id

    atom_str = sprint(show, res.atoms[1])
    atom = eval(Meta.parse(atom_str))
    @test atom == res.atoms[1]

    res_str = sprint(show, res)
    new_res = eval(Meta.parse(res_str))
    @test new_res == res
end
