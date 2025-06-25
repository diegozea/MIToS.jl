@testset "MMCIFDict conversions" begin
    cif_file = joinpath(DATA, "2vqc.cif")

    @testset "label = true" begin
        residues = read_file(cif_file, MMCIFFile)
        dict = BioStructures.MMCIFDict(residues; label = true)
        @test convert(BioStructures.MMCIFDict, residues) ==
              BioStructures.MMCIFDict(residues)
        @test convert(Vector{PDBResidue}, dict) == residues
        ref = BioStructures.MMCIFDict(cif_file)
        numeric_fields = Set([
            "_atom_site.Cartn_x",
            "_atom_site.Cartn_y",
            "_atom_site.Cartn_z",
            "_atom_site.occupancy",
            "_atom_site.B_iso_or_equiv",
        ])
        for field in keys(dict)
            @test haskey(ref, field)
            if field in numeric_fields
                @test all(
                    parse(Float64, d) == parse(Float64, r) for
                    (d, r) in zip(dict[field], ref[field])
                )
            else
                @test ref[field] == dict[field]
            end
        end
    end

    @testset "label = false" begin
        residues = read_file(cif_file, MMCIFFile, label = false)
        dict = BioStructures.MMCIFDict(residues; label = false)
        @test convert(Vector{PDBResidue}, dict) == residues
        @test convert(BioStructures.MMCIFDict, residues) ==
              BioStructures.MMCIFDict(residues)
        ref = BioStructures.MMCIFDict(cif_file)
        numeric_fields = Set([
            "_atom_site.Cartn_x",
            "_atom_site.Cartn_y",
            "_atom_site.Cartn_z",
            "_atom_site.occupancy",
            "_atom_site.B_iso_or_equiv",
        ])
        for field in keys(dict)
            @test haskey(ref, field)
            if field in numeric_fields
                @test all(
                    parse(Float64, d) == parse(Float64, r) for
                    (d, r) in zip(dict[field], ref[field])
                )
            else
                @test ref[field] == dict[field]
            end
        end
    end
end

@testset "MolecularStructure from residues" begin
    pdb_file = joinpath(DATA, "1SSX.pdb")
    residues = read_file(pdb_file, PDBFile)
    str_from_vec = BioStructures.MolecularStructure(residues)
    str_ref = BioStructures.read(pdb_file, BioStructures.PDBFormat)

    @test BioStructures.countmodels(str_from_vec) == BioStructures.countmodels(str_ref)
    @test BioStructures.countchains(str_from_vec) == BioStructures.countchains(str_ref)
    @test BioStructures.countresidues(str_from_vec) == BioStructures.countresidues(str_ref)
    @test BioStructures.countatoms(str_from_vec) == BioStructures.countatoms(str_ref)

    info(r) = (
        BioStructures.resname(r),
        BioStructures.resnumber(r),
        BioStructures.inscode(r),
        BioStructures.ishetero(r),
        length(r),
    )
    @test [info(r) for r in BioStructures.collectresidues(str_from_vec)] == [info(r) for r in BioStructures.collectresidues(str_ref)]

    atominfo(res) = [
        (
            BioStructures.atomname(a),
            BioStructures.altlocid(a) == '\0' ? ' ' : BioStructures.altlocid(a),
        ) for a in res
    ]
    @test [atominfo(r) for r in BioStructures.collectresidues(str_from_vec)] == [atominfo(r) for r in BioStructures.collectresidues(str_ref)]

    first_two = [
        (BioStructures.resnumber(r), BioStructures.inscode(r)) for
        r in BioStructures.collectresidues(str_from_vec)[1:2]
    ]
    @test first_two == [(15, 'A'), (15, 'B')]

    function altlocs_ca(str)
        for r in BioStructures.collectresidues(str)
            BioStructures.resnumber(r) == 90 || continue
            for at in r
                if BioStructures.atomname(at) == "CA"
                    return sort(BioStructures.altlocids(at))
                end
            end
        end
        Char[]
    end
    @test altlocs_ca(str_from_vec) == altlocs_ca(str_ref) == ['A', 'B']
end

@testset "Vector{PDBResidue} from MolecularStructure" begin
    pdb_file = joinpath(DATA, "1SSX.pdb")
    residues = read_file(pdb_file, PDBFile)
    str_ref = BioStructures.read(pdb_file, BioStructures.PDBFormat)
    str_from_vec = BioStructures.MolecularStructure(residues)

    conv_ref = convert(Vector{PDBResidue}, str_ref)
    conv_vec = convert(Vector{PDBResidue}, str_from_vec)

    @test length(conv_ref) == length(residues)
    @test length(conv_vec) == length(residues)

    recon = BioStructures.MolecularStructure(conv_ref)
    @test BioStructures.countmodels(recon) == BioStructures.countmodels(str_ref)
    @test BioStructures.countchains(recon) == BioStructures.countchains(str_ref)
    @test BioStructures.countresidues(recon) == BioStructures.countresidues(str_ref)
    @test BioStructures.countatoms(recon) == BioStructures.countatoms(str_ref)
end
