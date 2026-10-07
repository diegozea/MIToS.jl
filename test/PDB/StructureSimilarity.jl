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
            push!(ranges, offset:(offset+len-1))
            offset += len
        end
        perturbed = copy(coords)
        for r in ranges[(end-n+1):end]
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

@testset "Structure similarity contracts" begin
    function ca_residues(coords; group = "ATOM")
        [
            PDBResidue(
                PDBResidueIdentifier(string(i), string(i), "ALA", group, "1", "A"),
                [PDBAtom(Coordinates(coords[i, :]), "CA", "C", 1.0, "0", "", "")],
            ) for i in axes(coords, 1)
        ]
    end

    # These symmetric pairs have an identity least-squares fit and known distances.
    coords = [10.0 0 0; -10 0 0; 0 10 0; 0 -10 0; 0 0 10; 0 0 -10]
    reference = ca_residues(coords)
    model = ca_residues(coords .* [1.05 1.1 1.2])

    @testset "Definitions and numerical values" begin
        scores = gdt_per_cutoff(
            reference,
            model;
            cutoffs = (0.5, 1, 2, 4, 8),
            local_search = false,
        )
        @test scores[0.5] ≈ 100 / 3
        @test scores[1.0] ≈ 200 / 3
        @test scores[2.0] == scores[4.0] == scores[8.0] == 100.0
        @test gdt_ts(reference, model; local_search = false) ≈ 100 * 11 / 12
        @test gdt_ha(reference, model; local_search = false) ≈ 75.0
        @test tm_score(reference, model; local_search = false) ≈
              (1 / 2 + 1 / 5 + 1 / 17) / 3
        # Length changes both d0 and the denominator, not just an overall multiplier.
        d0 = 1.24 * cbrt(35) - 1.8
        expected = 2 * sum(1 / (1 + (d / d0)^2) for d in (0.5, 1.0, 2.0)) / 50
        @test tm_score(reference, model; Ltarget = 50, local_search = false) ≈ expected
        @test gdt_per_cutoff(reference, model; cutoffs = [2, 2], local_search = false) ==
              Dict(2.0 => 100.0)
        @test isempty(gdt_per_cutoff(reference, model; cutoffs = ()))
        for L in (1, 2, 3, 14, 15, 21, 22, 50)
            @test PDB._tm_d0(L) ≈ max(0.5, 1.24 * cbrt(Float64(L) - 15) - 1.8)
        end
        @test tm_score(reference[1:1], reference[1:1]; Ltarget = UInt(1)) ≈ 1.0
    end

    @testset "Reference normalization and supplied correspondence" begin
        for score in (gdt_ts, gdt_ha, tm_score), search in (false, true)
            perfect = score === tm_score ? 1.0 : 100.0
            @test score(reference, reference[1:3]; local_search = search) ≈ perfect / 2
            @test score(reference[1:3], reference; local_search = search) ≈ perfect
            @test score(reference, reference; Ltarget = 12, local_search = search) ≈
                  perfect / 2
            @test score(
                reference,
                reverse(reference);
                matches = ((i, 7 - i) for i = 1:6 if i != 2),
                local_search = search,
            ) ≈ perfect * 5 / 6
            @test score(reference, reference; matches = ((i, i) for i = 2:5)) ≈
                  perfect * 4 / 6
            @test score(reference, reference[[2, 1, 3, 4, 5, 6]]; local_search = search) <
                  perfect
        end
        @test gdt_ts(view(reference, 1:3), view(reference, 1:3)) ≈ 100.0
    end

    @testset "Iterative search recovers a core and retains all pair contributions" begin
        # An outlier translates the global centroid by 20/7 Å. The 4 Å GDT
        # refinement and the TM selection radius recover the six-residue core.
        target = ca_residues(vcat(coords, [0.0 0 0]))
        displaced = ca_residues(vcat(coords, [0.0 0 20]))
        for score in (gdt_ts, gdt_ha, tm_score)
            global_score = score(target, displaced; local_search = false)
            seed_score = score(target, displaced; window_sizes = (), max_iterations = 0)
            refined_score = score(target, displaced; window_sizes = ())
            @test global_score ≈ seed_score
            @test refined_score > seed_score
            @test score(target, displaced) >= refined_score
            @test score(target, displaced; max_iterations = 1) <= score(target, displaced)
        end
        @test gdt_ts(target, displaced; window_sizes = ()) ≈ 600 / 7
        @test gdt_ha(target, displaced; window_sizes = ()) ≈ 600 / 7
        expected_tm = (6 + 1 / (1 + (20 / 0.5)^2)) / 7
        @test tm_score(target, displaced; window_sizes = ()) ≈ expected_tm
        @test tm_score(target, displaced; window_sizes = ()) > 6 / 7
        # Selection must also make progress when fewer than three pairs are nearby.
        @test PDB._tm_subset([30.0, 10.0, 20.0, 40.0], 4.5) == [1, 2, 3]
        @test PDB._fragment_windows(6, [3, 4, 10]) ==
              [1:6, 1:3, 2:4, 3:5, 4:6, 1:4, 2:5, 3:6]
    end

    @testset "Cα selection and missing residues" begin
        missing_ca = deepcopy(reference)
        empty!(missing_ca[2].atoms)
        missing_ca[4].atoms = [PDBAtom(Coordinates(0, 0, 0), "N", "N", 1.0, "0", "", "")]
        for score in (gdt_ts, gdt_ha, tm_score)
            perfect = score === tm_score ? 1.0 : 100.0
            @test score(reference, missing_ca) ≈ perfect * 4 / 6
            @test score(missing_ca, reference) ≈ perfect * 4 / 6
            @test score(missing_ca, missing_ca) ≈ perfect * 4 / 6
        end
        modified = ca_residues(coords; group = "HETATM")
        @test tm_score(reference, modified) ≈ 1.0
        @test gdt_ts(reference, modified) ≈ 100.0
        alternative = deepcopy(reference)
        pushfirst!(
            alternative[1].atoms,
            PDBAtom(Coordinates(50, 50, 50), "CA", "C", 0.2, "0", "B", ""),
        )
        @test tm_score(reference, alternative) ≈ 1.0
        @test gdt_ha(reference, alternative) ≈ 100.0
        # Non-Cα atoms do not enter the fit and the caller's structures stay intact.
        push!(
            alternative[1].atoms,
            PDBAtom(Coordinates(NaN, 0, 0), "CB", "C", 1.0, "0", "", ""),
        )
        before_A, before_B = deepcopy(reference), deepcopy(alternative)
        @test tm_score(reference, alternative) ≈ 1.0
        @test gdt_ts(reference, alternative) ≈ 100.0
        @test isequal(reference, before_A)
        @test isequal(alternative, before_B)
    end

    @testset "Empty, short and degenerate structures" begin
        empty_structure = PDBResidue[]
        absent = [PDBResidue(reference[1].id, PDBAtom[])]
        for score in (gdt_ts, gdt_ha, tm_score)
            perfect = score === tm_score ? 1.0 : 100.0
            @test score(empty_structure, empty_structure) == 0.0
            @test score(reference, empty_structure) == 0.0
            @test score(empty_structure, reference) == 0.0
            @test score(absent, absent) == 0.0
            @test score(reference, reference; matches = ()) == 0.0
            for n = 1:3
                @test score(reference[1:n], reference[1:n]) ≈ perfect
            end
            line = ca_residues([0.0 0 0; 1 0 0; 2 0 0; 3 0 0])
            rotated = ca_residues([4.0 2 -1; 4 3 -1; 4 4 -1; 4 5 -1])
            @test score(line, rotated) ≈ perfect
            coincident = ca_residues(zeros(4, 3))
            @test score(coincident, coincident) ≈ perfect
        end
        @test gdt_per_cutoff(empty_structure, empty_structure) ==
              Dict(c => 0.0 for c in (1.0, 2.0, 4.0, 8.0))
    end

    @testset "Invalid arguments" begin
        for score in (gdt_ts, gdt_ha, tm_score)
            @test_throws BoundsError score(reference, model; matches = [(0, 1)])
            @test_throws BoundsError score(reference, model; matches = [(1, 7)])
            @test_throws ArgumentError score(reference, model; matches = [(1.0, 1)])
            for pairs in ([(1, 1), (1, 1)], [(1, 1), (1, 2)], [(1, 1), (2, 1)])
                @test_throws ArgumentError score(reference, model; matches = pairs)
            end
            for L in (-1, 0, 1)
                @test_throws ArgumentError score(reference, model; Ltarget = L)
            end
            @test_throws ArgumentError score(PDBResidue[], PDBResidue[]; Ltarget = -1)
            @test_throws ArgumentError score(reference, model; max_iterations = -1)
            for windows in ((0,), (-1,), (2,), (3.5,))
                @test_throws ArgumentError score(reference, model; window_sizes = windows)
            end
            for invalid in (NaN, Inf, -Inf), axis = 1:3
                invalid_coords = copy(coords)
                invalid_coords[1, axis] = invalid
                bad = ca_residues(invalid_coords)
                @test_throws ArgumentError score(reference, bad)
                @test_throws ArgumentError score(bad, reference)
            end
        end
        for cutoffs in ((-1,), (0,), (NaN,), (Inf,))
            @test_throws ArgumentError gdt_per_cutoff(reference, model; cutoffs = cutoffs)
        end
        @test_throws ArgumentError gdt_ts(reference, model; cutoffs = (1,))
        @test_throws ArgumentError gdt_ha(reference, model; cutoffs = (1,))
    end
end
