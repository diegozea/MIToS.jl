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
