@testset "PDBResidues" begin

    @test vec(Coordinates((1.0, 2.0, 3.0))) == [1.0, 2.0, 3.0]
    @test_throws ArgumentError Coordinates([1.0, 2.0])

    c1 = Coordinates(0.0, 0.0, 0.0)
    c2 = Coordinates(0.5, 0.5, sqrt(0.5))
    limit = 1.0
    @test contact(c1, c2, limit)
    c3 = Coordinates(0.5, 0.5, sqrt(0.5) + 0.001)
    @test !contact(c1, c3, limit)
    @test contact(c1, c1, limit)

    atom1 = PDBAtom(
        coordinates = c1,
        atom = "CA",
        element = "C",
        occupancy = 1.0,
        B = "0",
        alt_id = "",
        charge = "",
    )
    atom2 = PDBAtom(
        coordinates = c2,
        atom = "CA",
        element = "C",
        occupancy = 1.0,
        B = "0",
        alt_id = "",
        charge = "",
    )
    atom3 = PDBAtom(
        coordinates = c3,
        atom = "CA",
        element = "C",
        occupancy = 1.0,
        B = "0",
        alt_id = "",
        charge = "",
    )
    @test contact(atom1, atom2, limit)
    @test !contact(atom1, atom3, limit)
    @test contact(atom1, atom1, limit)

    res1 = PDBResidue(
        PDBResidueIdentifier("1", "1", "ALA", "ATOM", "1", "A"),
        [
            PDBAtom(
                coordinates = c1,
                atom = "CA",
                element = "C",
                occupancy = 1.0,
                B = "0",
                alt_id = "",
                charge = "",
            ),
            PDBAtom(
                coordinates = Coordinates(0.0, 3.0, 0.0),
                atom = "CB",
                element = "C",
                occupancy = 1.0,
                B = "0",
                alt_id = "",
                charge = "",
            ),
        ],
    )
    res2 = PDBResidue(
        PDBResidueIdentifier("2", "2", "ALA", "ATOM", "1", "A"),
        [
            PDBAtom(
                coordinates = Coordinates(1.0, 0.0, 0.0),
                atom = "CA",
                element = "C",
                occupancy = 1.0,
                B = "0",
                alt_id = "",
                charge = "",
            ),
            PDBAtom(
                coordinates = Coordinates(5.0, 0.0, 0.0),
                atom = "CB",
                element = "C",
                occupancy = 1.0,
                B = "0",
                alt_id = "",
                charge = "",
            ),
        ],
    )
    res3 = PDBResidue(
        PDBResidueIdentifier("3", "3", "ALA", "ATOM", "1", "A"),
        [
            PDBAtom(
                coordinates = Coordinates(10.0, 0.0, 0.0),
                atom = "CA",
                element = "C",
                occupancy = 1.0,
                B = "0",
                alt_id = "",
                charge = "",
            ),
            PDBAtom(
                coordinates = Coordinates(15.0, 0.0, 0.0),
                atom = "CB",
                element = "C",
                occupancy = 1.0,
                B = "0",
                alt_id = "",
                charge = "",
            ),
        ],
    )

    @test contact(res1, res2, 2.0)
    @test contact(res1, res2, 2.0, criteria = "CA")
    @test !contact(res1, res2, 2.0, criteria = "CB")
    @test !contact(res1, res3, 2.0)

    residues = [res1, res2, res3]
    expected = Bool[
        true true false;
        true true false;
        false false true
    ]
    @test contact(residues, 2.0) == expected

    @testset "print_file" begin
        res1 = read_file(joinpath(DATA, "short.pdb"), PDBFile)[1]

        io = IOBuffer()
        Utils.print_file(io, res1, PDBFile)
        str = String(take!(io))
        parsed = parse_file(IOBuffer(str), PDBFile)
        @test parsed == [res1]

        io = IOBuffer()
        Utils.print_file(io, res1, PDBFile, 10)
        lines = split(String(take!(io)), '\n'; keepempty = false)
        @test startswith(lines[1], "ATOM     10  N")
        @test startswith(lines[2], "ATOM     11  CA")
    end
end

@testset "cross PDBAtom" begin
    a = PDBAtom(
        coordinates = Coordinates(0.0, 1.0, 0.0),
        atom = "CA",
        element = "C",
        occupancy = 1.0,
        B = "0",
        alt_id = "",
        charge = "",
    )
    b = PDBAtom(
        coordinates = Coordinates(0.0, 0.0, 1.0),
        atom = "CB",
        element = "C",
        occupancy = 1.0,
        B = "0",
        alt_id = "",
        charge = "",
    )

    @test cross(a, b) == cross(a.coordinates, b.coordinates)
    @test cross(a, b) == Coordinates(1.0, 0.0, 0.0)
    @test cross(b, a) == Coordinates(-1.0, 0.0, 0.0)
end

@testset "change_b_factor" begin
    atom = PDBAtom(
        coordinates = Coordinates(0.0, 0.0, 0.0),
        atom = "CA",
        element = "C",
        occupancy = 1.0,
        B = "0.00",
        alt_id = "",
        charge = "",
    )
    newatom = change_b_factor(atom, 1.5)
    @test newatom.B == "1.50"
    @test atom.B == "0.00"  # original atom unchanged
    @test_throws ErrorException change_b_factor(atom, 123456.7)

    residue = PDBResidue(
        PDBResidueIdentifier("1", "1", "ALA", "ATOM", "1", "A"),
        [
            atom,
            PDBAtom(
                coordinates = Coordinates(1.0, 0.0, 0.0),
                atom = "CB",
                element = "C",
                occupancy = 1.0,
                B = "0.00",
                alt_id = "",
                charge = "",
            ),
        ],
    )
    newres = change_b_factor(residue, 2.0)
    @test all(a.B == "2.00" for a in newres.atoms)
    @test residue.atoms[1].B == "0.00"
    newres_ca = change_b_factor(residue, 3.0; atom = "CA")
    @test newres_ca.atoms[1].B == "3.00"
    @test newres_ca.atoms[2].B == "0.00"

    rescopy = deepcopy(residue)
    change_b_factor!(rescopy, 4.0)
    @test all(a.B == "4.00" for a in rescopy.atoms)
    @test residue.atoms[1].B == "0.00"

    rescopy2 = deepcopy(residue)
    change_b_factor!(rescopy2, 5.0; atom = "CA")
    @test rescopy2.atoms[1].B == "5.00"
    @test rescopy2.atoms[2].B == "0.00"
end
