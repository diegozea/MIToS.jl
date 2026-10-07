@testset "AlphaFold mmCIF templates" begin
    residue(number; name = "ALA", chain = "AA", model = "7", group = "ATOM", pdbe = "") =
        PDBResidue(
            PDBResidueIdentifier(pdbe, number, name, group, model, chain),
            [PDBAtom(Coordinates(1.0, 2.0, 3.0), "CA", "C", 1.0, "12.3", "", "")],
        )
    template(residues; kwargs...) =
        alphafold_mmcifdict(residues; entry_id = "local_template", kwargs...)

    # Nonmonotonic author numbers, reverse insertion order, interleaved chains and
    # inconsistent stored PDBe numbers must not dictate the new polymer sequence.
    residues = [
        residue("10B"; pdbe = "99"),
        residue("-2"; name = "SER", chain = "Z"),
        residue("10A"; name = "GLY", pdbe = "99"),
        residue("0"; name = "THR", chain = "Z"),
        residue("10"; name = "VAL"),
    ]
    push!(residues[1].atoms, PDBAtom(Coordinates(2, 3, 4), "N", "N", 0.5, "8", "B", "1"))
    original = deepcopy(residues)
    mmcif = template(view(residues, :))

    @testset "Identifiers and category relationships" begin
        @test mmcif["data_"] == mmcif["_entry.id"] == ["local_template"]
        @test mmcif["_entity.id"] == ["1", "2"]
        @test mmcif["_entity.type"] == ["polymer", "polymer"]
        @test mmcif["_struct_asym.id"] == ["AA", "Z"]
        @test mmcif["_struct_asym.entity_id"] == ["1", "2"]
        @test mmcif["_entity_poly_seq.entity_id"] == ["1", "1", "1", "2", "2"]
        @test mmcif["_entity_poly_seq.num"] == ["1", "2", "3", "1", "2"]
        @test mmcif["_entity_poly_seq.mon_id"] == ["ALA", "GLY", "VAL", "SER", "THR"]
        @test mmcif["_atom_site.label_seq_id"] == ["1", "1", "2", "3", "1", "2"]
        @test mmcif["_atom_site.auth_seq_id"] == ["10", "10", "10", "10", "-2", "0"]
        @test mmcif["_atom_site.pdbx_PDB_ins_code"] == ["B", "B", "A", "?", "?", "?"]
        @test mmcif["_atom_site.label_entity_id"] == ["1", "1", "1", "1", "2", "2"]
        @test mmcif["_atom_site.id"] == string.(1:6)
        @test all(==("1"), mmcif["_atom_site.pdbx_PDB_model_num"])
        for field in ("asym_id", "comp_id", "atom_id")
            @test mmcif["_atom_site.label_$field"] == mmcif["_atom_site.auth_$field"]
        end
        # Every atom refers to the matching entity, sequence position and monomer.
        seq = Dict(
            zip(
                zip(mmcif["_entity_poly_seq.entity_id"], mmcif["_entity_poly_seq.num"]),
                mmcif["_entity_poly_seq.mon_id"],
            ),
        )
        for (entity, number, name) in zip(
            mmcif["_atom_site.label_entity_id"],
            mmcif["_atom_site.label_seq_id"],
            mmcif["_atom_site.label_comp_id"],
        )
            @test seq[(entity, number)] == name
        end
        @test all(
            length(value) == 6 for (key, value) in mmcif if startswith(key, "_atom_site.")
        )
        @test residues == original
    end

    @testset "Serialization and round trip" begin
        mktemp() do path, io
            close(io)
            BioStructures.writemmcif(path, mmcif)
            @test BioStructures.MMCIFDict(path) == mmcif
            parsed = read_file(path, MMCIFFile)
            expected = original[[1, 3, 5, 2, 4]]
            @test [r.id.number for r in parsed] == [r.id.number for r in expected]
            @test [r.id.chain for r in parsed] == [r.id.chain for r in expected]
            @test [r.atoms for r in parsed] == [r.atoms for r in expected]
            @test read_file(path, MMCIFFile; label = false) == parsed
        end
        renumbered = template(residues; renumber_auth = true)
        @test renumbered["_atom_site.auth_seq_id"] == renumbered["_atom_site.label_seq_id"]
        @test all(==("?"), renumbered["_atom_site.pdbx_PDB_ins_code"])
        @test residues == original
    end

    @testset "Model selection" begin
        models = vcat(residues, [residue("1"; model = "2", name = "MET")])
        @test_throws ArgumentError template(models)
        @test_throws ArgumentError template(models; model = "3")
        @test template(models; model = "7") == mmcif
        selected = template(models; model = "2")
        @test selected["_entity_poly_seq.mon_id"] == ["MET"]
        @test selected["_entity_poly_seq.num"] == ["1"]
        @test selected["_atom_site.pdbx_PDB_model_num"] == ["1"]
    end

    @testset "Scientific metadata is explicit" begin
        types = Dict(zip(mmcif["_chem_comp.id"], mmcif["_chem_comp.type"]))
        @test types["ALA"] == "L-peptide linking"
        @test types["GLY"] == "peptide linking"
        @test mmcif["_exptl.method"] == ["?"]
        @test !haskey(mmcif, "_pdbx_audit_revision_history.revision_date")
        @test !any(
            occursin("resolution", key) || startswith(key, "_refine.") for
            key in keys(mmcif)
        )
        known = template(
            residues;
            exptl_method = "SOLUTION NMR",
            release_date = PDB.Date(2020, 2, 3),
        )
        @test known["_exptl.method"] == ["SOLUTION NMR"]
        @test known["_pdbx_audit_revision_history.revision_date"] == ["2020-02-03"]

        for name in ("MSE", "UNK", "DAL", "HOH", "DA", "ATP", "ZZZ")
            @test_throws ArgumentError template([
                residue("1"; name = name, group = "HETATM"),
            ])
        end
        modified = template(
            [residue("1"; name = "MSE", group = "HETATM"), residue("2"; name = "DAL")];
            chem_comp_types = Dict(
                "MSE" => "L-peptide linking",
                "DAL" => "D-peptide linking",
            ),
        )
        @test modified["_chem_comp.type"] == ["L-peptide linking", "D-peptide linking"]
        @test modified["_entity_poly_seq.mon_id"] == ["MSE", "DAL"]
        @test modified["_atom_site.group_PDB"] == ["HETATM", "ATOM"]
        for value in ("non-polymer", "not a peptide", "", "?", 42)
            @test_throws ArgumentError template(
                residues;
                chem_comp_types = Dict("ALA" => value),
            )
        end
        missing_element = [residue("1")]
        missing_element[1].atoms[1] =
            PDBAtom(Coordinates(0, 0, 0), "CA", "", 1.0, "0", "", "")
        @test template(missing_element)["_atom_site.type_symbol"] == ["?"]
    end

    @testset "Reject ambiguous input" begin
        @test_throws ArgumentError template(PDBResidue[])
        for chain in ("", " ", "?", ".", "A B")
            @test_throws ArgumentError template([residue("1"; chain = chain)])
        end
        for number in (
            "",
            ".",
            "?",
            "ABC",
            "10AB",
            "10.1",
            "A10",
            "10 A",
            "10A\n",
            "999999999999999999999999",
        )
            @test_throws ArgumentError template([residue(number)])
        end
        @test template([residue("+001A")])["_atom_site.auth_seq_id"] == ["1"]
        for duplicate in (residue("10B"), residue("010B"; name = "SER", group = "HETATM"))
            @test_throws ArgumentError template(vcat(residues, [duplicate]))
            @test_throws ArgumentError template(
                vcat(residues, [duplicate]);
                renumber_auth = true,
            )
        end
        empty_residue = residue("1")
        empty!(empty_residue.atoms)
        @test_throws ArgumentError template([empty_residue])
        @test_throws ArgumentError template([residue("1"; group = "OTHER")])
        for name in ("", ".", "?", "A B")
            @test_throws ArgumentError template(
                [residue("1"; name = name)];
                chem_comp_types = Dict(name => "peptide linking"),
            )
        end
        for entry in ("", "?", ".", "two words", "bad\nentry", "entry\n")
            @test_throws ArgumentError alphafold_mmcifdict(residues; entry_id = entry)
        end
        @test_throws ArgumentError template(residues; exptl_method = " ")
    end

    @testset "Generic export remains atom_site only" begin
        for label in (true, false)
            generic = BioStructures.MMCIFDict(original; label = label)
            @test all(startswith(key, "_atom_site.") for key in keys(generic))
            @test generic["_atom_site.label_seq_id"] == ["99", "99", ".", "99", ".", "."]
            @test all(==("7"), generic["_atom_site.pdbx_PDB_model_num"])
            @test haskey(generic, "_atom_site.label_asym_id") == label
            @test haskey(generic, "_atom_site.auth_asym_id") == !label
            io = IOBuffer()
            print_file(io, original, MMCIFFile; label = label)
            @test parse_file(IOBuffer(take!(io)), MMCIFFile; label = label) == original
        end
    end

    @testset "PDB input regression" begin
        pdb = read_file(joinpath(DATA, "1SSX.pdb"), PDBFile; chain = "A", group = "ATOM")
        @test all(isempty(r.id.PDBe_number) for r in pdb)
        cif = template(pdb)
        @test cif["_entity_poly_seq.num"] == string.(1:length(pdb))
        @test cif["_entity_poly_seq.mon_id"] == [r.id.name for r in pdb]
        @test cif["_atom_site.pdbx_PDB_ins_code"][1] == "A"
        nmr = read_file(joinpath(DATA, "1AS5.pdb"), PDBFile)
        # The terminal NH2 cap is not an amino acid; select the peptide residues.
        filter!(r -> r.id.name != "NH2", nmr)
        @test_throws ArgumentError template(nmr)
        selected = template(
            nmr;
            model = "14",
            chem_comp_types = Dict("HYP" => "L-peptide linking"),
        )
        @test length(selected["_entity_poly_seq.num"]) == 24
        @test all(==("1"), selected["_atom_site.pdbx_PDB_model_num"])
    end
end
