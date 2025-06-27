@testset "MMCIF parse_file filters" begin

    occ_res = PDBResidue(
        PDBResidueIdentifier("1", "1", "ALA", "ATOM", "1", "A"),
        [
            PDBAtom(
                coordinates = Coordinates(0.0, 0.0, 0.0),
                atom = "CA",
                element = "C",
                occupancy = 0.4,
                B = "0",
                alt_id = "A",
                charge = "",
            ),
            PDBAtom(
                coordinates = Coordinates(1.0, 0.0, 0.0),
                atom = "CA",
                element = "C",
                occupancy = 0.6,
                B = "0",
                alt_id = "B",
                charge = "",
            ),
        ],
    )

    heavy_res = PDBResidue(
        PDBResidueIdentifier("2", "2", "GLY", "ATOM", "1", "A"),
        [
            PDBAtom(
                coordinates = Coordinates(0.0, 0.0, 0.0),
                atom = "CA",
                element = "C",
                occupancy = 1.0,
                B = "0",
                alt_id = "",
                charge = "",
            ),
            PDBAtom(
                coordinates = Coordinates(1.0, 0.0, 0.0),
                atom = "H1",
                element = "H",
                occupancy = 1.0,
                B = "0",
                alt_id = "",
                charge = "",
            ),
        ],
    )

    residues = [occ_res, heavy_res]
    io = IOBuffer()
    Utils.print_file(io, residues, MMCIFFile; label = true)
    cif_str = String(take!(io))

    # Without filters
    parsed = parse_file(IOBuffer(cif_str), MMCIFFile)
    @test length(parsed[1].atoms) == 2
    @test length(parsed[2].atoms) == 2

    # Occupancy filter
    parsed_occ = parse_file(IOBuffer(cif_str), MMCIFFile; occupancyfilter = true)
    @test length(parsed_occ[1].atoms) == 1
    @test parsed_occ[1].atoms[1].occupancy == 0.6
    @test length(parsed_occ[2].atoms) == 2

    # Only heavy atoms
    parsed_heavy = parse_file(IOBuffer(cif_str), MMCIFFile; onlyheavy = true)
    @test length(parsed_heavy[1].atoms) == 2
    @test length(parsed_heavy[2].atoms) == 1
    @test parsed_heavy[2].atoms[1].element == "C"
end
