@testset "Copy and deepcopy" begin
    pdbfile = joinpath(DATA, "short.pdb")
    residue = read_file(pdbfile, PDBFile)[1]

    # immutable objects return the same instance
    id_copy = copy(residue.id)
    id_deepcopy = deepcopy(residue.id)
    @test id_copy === residue.id
    @test id_deepcopy === residue.id

    coords = residue.atoms[1].coordinates
    coords_copy = copy(coords)
    coords_deepcopy = deepcopy(coords)
    @test coords_copy === coords
    @test coords_deepcopy === coords

    atom = residue.atoms[1]
    atom_copy = copy(atom)
    atom_deepcopy = deepcopy(atom)
    @test atom_copy === atom
    @test atom_deepcopy === atom

    res_copy = copy(residue)
    res_deepcopy = deepcopy(residue)
    @test res_copy == residue
    @test res_deepcopy == residue
    @test res_copy !== residue
    @test res_deepcopy !== residue

    # modifications in the copy don't change the original
    n = length(residue.atoms)
    push!(res_copy.atoms, copy(res_copy.atoms[1]))
    @test length(res_copy.atoms) == n + 1
    @test length(residue.atoms) == n
    @test length(res_deepcopy.atoms) == n

    # modifications in the original don't change the copy
    pop!(residue.atoms)
    @test length(residue.atoms) == n - 1
    @test length(res_copy.atoms) == n + 1
    @test length(res_deepcopy.atoms) == n
end
