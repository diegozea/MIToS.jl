@testset "Structure similarity (GDT/TM-score)" begin
    structure = read_file(joinpath(DATA, "1CBN.pdb"), PDBFile, group = "ATOM", model = "1")

    @test gdt_ts(structure, structure) ≈ 100.0 atol = 1.0e-8
    @test gdt_ha(structure, structure) ≈ 100.0 atol = 1.0e-8
    @test tm_score(structure, structure) ≈ 1.0 atol = 1.0e-8

    rot_z = [0.0 -1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 1.0]
    trans = [3.0, -2.0, 1.5]

    function transform_residues(residues, rotation, translation)
        coords = coordinatesmatrix(residues)
        transformed = coords * rotation .+ transpose(translation)
        return change_coordinates(residues, transformed)
    end

    transformed = transform_residues(structure, rot_z, trans)

    @test gdt_ts(structure, transformed) ≈ 100.0 atol = 1.0e-6
    @test gdt_ha(structure, transformed) ≈ 100.0 atol = 1.0e-6
    @test tm_score(structure, transformed) ≈ 1.0 atol = 1.0e-6
    @test gdt_ts(structure, transformed; local_search = false) ≈ 100.0 atol = 1.0e-6

    function perturb_last_residues(residues, n::Int, shift)
        coords = coordinatesmatrix(residues)
        offset = 1
        ranges = Vector{UnitRange{Int}}()
        for res in residues
            len = length(res)
            push!(ranges, offset:(offset + len - 1))
            offset += len
        end
        perturbed = copy(coords)
        for r in ranges[(end - n + 1):end]
            perturbed[r, :] .+= reshape(shift, 1, :)
        end
        return change_coordinates(residues, perturbed)
    end

    shifted = perturb_last_residues(structure, 5, [5.0, 0.0, 0.0])

    noisy_gdt = gdt_ts(structure, shifted)
    noisy_ha = gdt_ha(structure, shifted)
    noisy_tm = tm_score(structure, shifted)

    @test noisy_gdt < 100.0
    @test noisy_ha < noisy_gdt
    @test 0.0 < noisy_tm < 1.0
end
