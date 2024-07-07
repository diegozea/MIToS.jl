@testset "Parse PDB and PDBML" begin

    txt(code) = joinpath(DATA, string(uppercase(code), ".pdb"))
    xml(code) = joinpath(DATA, string(uppercase(code), ".xml"))

    @testset "2VQC: Missings" begin

        code = "2VQC"
        pdb = read_file(txt(code), PDBFile)
        pdbml = read_file(xml(code), PDBML)

        @test findfirst(x -> x.id.number == "4", pdb) ==
              findfirst(x -> x.id.number == "4", pdbml)
        @test findfirst(x -> x.id.number == "73", pdb) ==
              findfirst(x -> x.id.number == "73", pdbml)
    end

    @testset "Test downloadpdb" begin

        code = "2VQC"
        pdb = read_file(txt(code), PDBFile) # reference
        # PDBML
        filename = downloadpdb(code, format = PDBML)
        try
            pdbml = read_file(filename, PDBML)

            @test findfirst(x -> x.id.number == "4", pdb) ==
                  findfirst(x -> x.id.number == "4", pdbml)
            @test findfirst(x -> x.id.number == "73", pdb) ==
                  findfirst(x -> x.id.number == "73", pdbml)
        finally
            rm(filename)
        end
        # mmCIF (the default format)
        filename = downloadpdb(code)
        try
            mmcif = read_file(filename, MMCIFFile)

            @test findfirst(x -> x.id.number == "4", pdb) ==
                  findfirst(x -> x.id.number == "4", mmcif)
            @test findfirst(x -> x.id.number == "73", pdb) ==
                  findfirst(x -> x.id.number == "73", mmcif)
        finally
            rm(filename)
        end
        # PDB
        filename = downloadpdb(
            code,
            headers = Dict(
                "User-Agent" => "Mozilla/5.0 (compatible; MSIE 7.01; Windows NT 5.0)",
            ),
            format = PDBFile,
        )
        try
            d_pdb = read_file(filename, PDBFile)

            @test findfirst(x -> x.id.number == "4", pdb) ==
                  findfirst(x -> x.id.number == "4", d_pdb)
            @test findfirst(x -> x.id.number == "73", pdb) ==
                  findfirst(x -> x.id.number == "73", d_pdb)
        finally
            rm(filename)
        end
    end

    @testset "1H4A: Chain A (auth) == Chain X (label)" begin

        file = string(xml("1H4A"), ".gz")
        auth = read_file(file, PDBML, label = false)
        label = read_file(file, PDBML)

        @test unique([res.id.chain for res in auth]) == ["X"]
        @test unique([res.id.chain for res in label]) == ["A", "B"]
    end

    @testset "1SSX: Residues with insert codes, 15A 15B" begin

        code = "1SSX"
        pdb = read_file(txt(code), PDBFile)
        pdbml = read_file(xml(code), PDBML)

        for residue_list in [pdb, pdbml]

            @test findall(res -> res.id.number == "15A", residue_list) == [1]
            @test findall(res -> isresidue(res, residue = "15A"), residue_list) == [1]

            @test findall(res -> res.id.number == "15B", residue_list) == [2]
            @test findall(res -> isresidue(res, residue = "15B"), residue_list) == [2]
        end

        @testset "Occupancy != 1.0" begin

            @test sum(
                map(
                    a -> (a.atom == "HH22" ? a.occupancy : 0.0),
                    filter(r -> r.id.number == "141", pdbml)[1].atoms,
                ),
            ) == 1.0
            @test sum([
                a.occupancy for a in select_atoms(
                    pdbml,
                    model = "1",
                    chain = "A",
                    group = All,
                    residue = "141",
                    atom = "HH22",
                )
            ]) == 1.0
        end

        @testset "Best occupancy" begin

            atoms_141 = select_atoms(
                pdbml,
                model = "1",
                chain = "A",
                group = All,
                residue = "141",
                atom = "HH22",
            )
            resid_141 = select_residues(
                pdbml,
                model = "1",
                chain = "A",
                group = All,
                residue = "141",
            )

            @test bestoccupancy(atoms_141)[1].occupancy == 0.75
            @test bestoccupancy(reverse(atoms_141))[1].occupancy == 0.75
            @test bestoccupancy(PDBAtom[atoms_141[2]])[1].occupancy == 0.25

            @test length(resid_141[1]) == 48
            @test selectbestoccupancy(resid_141[1], collect(1:48)) == 1
            @test selectbestoccupancy(resid_141[1], [1, 2]) == 1

            @test_throws AssertionError selectbestoccupancy(resid_141[1], Int[])
            @test_throws AssertionError selectbestoccupancy(resid_141[1], collect(1:100))
        end

        @testset "select_atom with All" begin

            # ATOM      2  CA  ALA A  15A     22.554  11.619   6.400  1.00  6.14           C
            @test select_atoms(
                pdb,
                model = "1",
                chain = "A",
                group = "ATOM",
                residue = All,
                atom = r"C.+",
            )[1].atom == "CA"
        end
    end

    @testset "read_file with occupancyfilter=true" begin
        # `read_file` only atoms with the best occupancy

        code = "1SSX"
        pdb = read_file(txt(code), PDBFile, occupancyfilter = true)
        pdbml = read_file(xml(code), PDBML, occupancyfilter = true)

        res_pdb =
            select_residues(pdb, model = "1", chain = "A", group = All, residue = "141")
        res_pdbml =
            select_residues(pdbml, model = "1", chain = "A", group = All, residue = "141")

        atm_pdbml = select_atoms(
            pdbml,
            model = "1",
            chain = "A",
            group = All,
            residue = "141",
            atom = "HH22",
        )

        @test length(atm_pdbml) == 1
        @test atm_pdbml[1].occupancy == 0.75
        @test length(
            select_atoms(
                pdb,
                model = "1",
                chain = "A",
                group = All,
                residue = "141",
                atom = "HH22",
            ),
        ) == 1

        @test length(res_pdb[1]) == 24
        @test length(res_pdbml[1]) == 24
    end

    @testset "CBN: Identical PDBe ResNum for Residue 22" begin
        # <residue dbSource="PDBe" dbCoordSys="PDBe" dbResNum="22" dbResName="SER">
        #   <crossRefDb dbSource="PDB" dbCoordSys="PDBresnum" dbAccessionId="1cbn" dbResNum="22" dbResName="SER" dbChainId="A"/>
        #   <crossRefDb dbSource="UniProt" dbCoordSys="UniProt" dbAccessionId="P01542" dbResNum="22" dbResName="P"/>
        #   ....
        # <residue dbSource="PDBe" dbCoordSys="PDBe" dbResNum="22" dbResName="PRO">
        #   <crossRefDb dbSource="PDB" dbCoordSys="PDBresnum" dbAccessionId="1cbn" dbResNum="22" dbResName="PRO" dbChainId="A"/>
        #   <crossRefDb dbSource="UniProt" dbCoordSys="UniProt" dbAccessionId="P01542" dbResNum="22" dbResName="P"/>
        #   ...

        code = "1CBN"
        pdb = read_file(txt(code), PDBFile)
        pdbml = read_file(xml(code), PDBML)

        @test length(
            select_residues(pdb, model = "1", chain = "A", group = All, residue = "22"),
        ) == 2
        @test length(
            select_residues(pdbml, model = "1", chain = "A", group = All, residue = "22"),
        ) == 2

        @test [
            r.id.name for
            r in select_residues(pdb, model = "1", chain = "A", group = All, residue = "22")
        ] == ["SER", "PRO"]
        @test [
            r.id.name for r in
            select_residues(pdbml, model = "1", chain = "A", group = All, residue = "22")
        ] == ["SER", "PRO"]
    end

    @testset "1AS5: NMR" begin

        code = "1AS5"
        pdb = read_file(txt(code), PDBFile)
        pdbml = read_file(xml(code), PDBML)

        @test length(select_residues(pdbml, model = "1", chain = "A")) == 25
        @test length(select_residues(pdbml, model = "14", chain = "A")) == 25

        @test length(select_residues(pdbml, model = All, chain = "A")) == 25 * 14
    end

    @testset "1DPO: Inserted residues lack insertion letters" begin
        # Single unnamed chain in 1DPO contains insertions at postions 184 (Gly, Phe),
        # 188 (Gly, Lys), and 221 (Ala, Leu) but no insertion letters.

        code = "1DPO"
        pdb = read_file(txt(code), PDBFile)
        pdbml = read_file(xml(code), PDBML)

        # Single "A" chain for PDB (auth_asym_id in PDBML)
        @test unique([r.id.chain for r in pdb]) == ["A"]
        # But 'A':'H' chains for PDBML (label_asym_id)
        @test unique([r.id.chain for r in pdbml]) == [string(chain) for chain = 'A':'H']

        @test [
            r.id.name for
            r in select_residues(pdb, model = "1", chain = "A", residue = r"^184[A-Z]?$")
        ] == ["GLY", "PHE"]
        @test [
            r.id.name for
            r in select_residues(pdbml, model = "1", chain = "A", residue = r"^184[A-Z]?$")
        ] == ["GLY", "PHE"]

        @test [
            r.id.name for
            r in select_residues(pdb, model = "1", chain = "A", residue = r"^188[A-Z]?$")
        ] == ["GLY", "LYS"]
        @test [
            r.id.name for
            r in select_residues(pdbml, model = "1", chain = "A", residue = r"^188[A-Z]*")
        ] == ["GLY", "LYS"]

        @test [
            r.id.name for
            r in select_residues(pdb, model = "1", chain = "A", residue = r"^221[A-Z]?$")
        ] == ["ALA", "LEU"]
        @test [
            r.id.name for
            r in select_residues(pdbml, model = "1", chain = "A", residue = r"^221[A-Z]?$")
        ] == ["ALA", "LEU"]
    end

    @testset "1IGY: Insertions" begin
        # Insertions have more than one copy of the same amino acid in a single insertion block.
        # For example, chain B in 1IGY contains a block of four residues inserted at sequence position 82.
        # The block contains Leu-Ser-Ser-Leu.

        code = "1IGY"
        pdb = read_file(txt(code), PDBFile)
        pdbml = read_file(xml(code), PDBML)

        @test [
            r.id.name for r in select_residues(
                pdb,
                model = "1",
                chain = "B",
                group = All,
                residue = r"^82[A-Z]?$",
            )
        ] == ["LEU", "SER", "SER", "LEU"]
        @test [
            r.id.name for r in select_residues(
                pdbml,
                model = "1",
                chain = "B",
                group = All,
                residue = r"^82[A-Z]?$",
            )
        ] == ["LEU", "SER", "SER", "LEU"]

        @test sum(
            [
                r.id.group for r in
                select_residues(pdb, model = "1", chain = "D", group = All, residue = All)
            ] .== "HETATM",
        ) == length(
            select_residues(pdb, model = "1", chain = "D", group = "HETATM", residue = All),
        )
        @test sum(
            [
                r.id.group for r in
                select_residues(pdb, model = "1", chain = "D", group = All, residue = All)
            ] .== "ATOM",
        ) == length(
            select_residues(pdb, model = "1", chain = "D", group = "ATOM", residue = All),
        )
    end

    @testset "1HAG" begin
        # Chain E begins with 1H, 1G, 1F, ... 1A, then 1 (in reverse alphabetic order)

        code = "1HAG"
        pdb = read_file(txt(code), PDBFile)
        pdbml = read_file(xml(code), PDBML)

        @test unique([res.id.chain for res in pdb]) == ["E", "I"]
        @test unique([res.id.chain for res in pdbml]) == ["A", "B", "C", "D", "E"]

        # The chain E of PDB is the chain A of PDBML
        @test [
            r.id.number for r in select_residues(
                pdb,
                model = "1",
                chain = "E",
                group = All,
                residue = r"^1[A-Z]?$",
            )
        ] == [string(1, code) for code in vcat(collect('H':-1:'A'), "")]
        @test [
            r.id.number for r in select_residues(
                pdbml,
                model = "1",
                chain = "A",
                group = All,
                residue = r"^1[A-Z]?$",
            )
        ] == [string(1, code) for code in vcat(collect('H':-1:'A'), "")]
    end

    @testset "1NSA" begin
        # Contains a single (unnamed) protein chain with sequence 7A-95A that continues 4-308.

        code = "1NSA"
        pdb = read_file(txt(code), PDBFile)
        pdbml = read_file(xml(code), PDBML)

        # Single "A" chain for PDB (auth_asym_id in PDBML)
        @test unique([r.id.chain for r in pdb]) == ["A"]
        # But 'A':'F' chains for PDBML (label_asym_id)
        @test unique([r.id.chain for r in pdbml]) == [string(chain) for chain = 'A':'F']

        ind = findall(r -> r.id.number == "95A", pdbml)[1]
        @test pdbml[ind+1].id.number == "4"

        ind = findall(r -> r.id.number == "95A", pdb)[1]
        @test pdb[ind+1].id.number == "4"
    end

    @testset "1IAO" begin
        # Contains in chain B (in this order) 1S, 2S, 323P-334P, 6-94, 94A, 95-188, 1T, 2T

        code = "1IAO"
        pdb = read_file(txt(code), PDBFile)
        pdbml = read_file(xml(code), PDBML)
        pdb_B = select_residues(pdb, model = "1", chain = "B")
        pdbml_B = select_residues(pdbml, model = "1", chain = "B")

        for B in [pdb_B, pdbml_B]
            @test B[findall(r -> r.id.number == "2S", B)[1]+1].id.number == "323P"
            @test B[findall(r -> r.id.number == "334P", B)[1]+1].id.number == "6"
            @test B[findall(r -> r.id.number == "94", B)[1]+1].id.number == "94A"
            @test B[findall(r -> r.id.number == "94A", B)[1]+1].id.number == "95"
            @test B[findall(r -> r.id.number == "188", B)[1]+1].id.number == "1T"
            @test B[findall(r -> r.id.number == "1T", B)[1]+1].id.number == "2T"
        end
    end

    @testset "3BTT" begin
        # ASN 115 from chain E has no CA

        code = "3BTT"
        pdb = read_file(txt(code), PDBFile)
        pdbml = read_file(xml(code), PDBML)

        res = pdb[97]
        res_ml = pdbml[97]

        res.id.name == "ASN"
        res_ml.id.name == "ASN"

        @test getCA(res) === missing
        @test getCA(res_ml) === missing
    end

    @testset "Foldseek" begin
        # PDB files generated by Foldseek have only 66 columns; element identifiers 
        # are missing (represented as empty strings). Those files only contain CA atoms.
        pdb_file = joinpath(DATA, "foldseek_example.pdb")
        pdb = read_file(pdb_file, PDBFile)

        # The example file has only 8 residues and one chain (A)
        for (i, res) in zip(1:8, pdb)
            @test res.id.group == "ATOM"
            @test res.id.number == string(i)
            @test res.id.chain == "A"
            @test length(res.atoms) == 1
            atom = res.atoms[1]
            @test atom.atom == "CA"
            @test isempty(atom.element) # ""
            # All residues in the example have 1.0 occupancy and 0.00 temperature factor
            @test atom.occupancy == 1.0
            @test atom.B == "0.00"
        end
    end
end

@testset "RESTful PDB Interface" begin

    @testset "Percent-encoding" begin
        @test PDB._escape_url_query("name=John Snow") == "name%3DJohn%20Snow"

        unreserved = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-_.~"
        @test PDB._escape_url_query(unreserved) == unreserved

        @test PDB._escape_url_query("£") == "%C2%A3"
        @test PDB._escape_url_query("€") == "%E2%82%AC"
    end

    @test getpdbdescription("4HHB")["rcsb_entry_info"]["resolution_combined"][1] == 1.74
    @test getpdbdescription("104D")["rcsb_entry_info"]["resolution_combined"] === nothing

    mktemp() do path, io
        filename = downloadpdbheader("4HHB", filename = path)
        @test isfile(filename)
        @test filename == path
        file_content = read(filename, String)
        @test occursin("resolution_combined", file_content)
        @test occursin("1.74", file_content)
    end
end

@testset "Write PDB files" begin

    txt(code) = joinpath(DATA, string(uppercase(code), ".pdb"))

    @testset "2VQC" begin
        code = "2VQC"
        io = IOBuffer()
        pdb = read_file(txt(code), PDBFile)
        print_file(io, pdb, PDBFile)
        printed = split(String(take!(io)), '\n')

        @test length(printed) == 609 # Only ATOM, HETATM & END + 1 because the trailing \n
        @test printed[1] ==
              "ATOM      1  N   THR A   4       2.431  19.617   6.520  1.00 24.37           N  "
        @test printed[607] ==
              "HETATM  607  O   HOH A2025      13.807  38.993   2.453  1.00 33.00           O  "
    end

    @testset "NMR" begin
        code = "1AS5"
        io = IOBuffer()

        pdb = read_file(txt(code), PDBFile)
        print_file(io, pdb, PDBFile)
        printed = split(String(take!(io)), '\n')

        @test sum(map(x -> startswith(x, "MODEL "), printed)) == 14 # 14 models
        @test sum(map(x -> x == "ENDMDL", printed)) == 14 # 14 models
    end

    @testset "2 Chains" begin
        code = "1IAO"
        io = IOBuffer()

        pdb = read_file(txt(code), PDBFile)
        print_file(io, pdb, PDBFile)
        printed = split(String(take!(io)), '\n')

        # MIToS only prints TER for the ATOM group if the chain changes.
        # Some modified residues are annotated as HETATM in the middle of the ATOM chain:
        # TER can not be printed from ATOM to HETATM if the chain doesn’t change.

        # Only prints TER between chain A and B
        @test sum(map(x -> startswith(x, "TER "), printed)) == 1

        @test filter!(s -> occursin(r"TER ", s), printed)[1] ==
              "TER    1418      TRP A 178 "
    end

    @testset "read_file/write consistency" begin
        io = IOBuffer()
        for code in ["2VQC", "1IAO", "1NSA", "1HAG", "1IGY", "1DPO", "1AS5", "1CBN", "1SSX"]

            readed = read_file(txt(code), PDBFile)
            print_file(io, readed, PDBFile)
            readed_writed_readed = parse_file(String(take!(io)), PDBFile)

            @test readed_writed_readed == readed
        end
    end
end
