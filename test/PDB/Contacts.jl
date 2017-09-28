@testset "Contacts" begin

    txt(code) = joinpath(pwd(), "data", string(uppercase(code), ".pdb"))

    @testset "Piccolo" begin
    # Using data from http://www-cryst.bioc.cam.ac.uk/~richard/piccolo/piccolo.php?PDB=1IGY (28/Sep/2015)

        code = "1IGY"
        pdb = read(txt(code), PDBFile)

        # Modify the next line if ligands are added to AtomsData.jl
        @test sum(check_atoms_for_interactions(r) for r in pdb) ==
              sum(r.id.group == "ATOM" for r in pdb)

        @testset "findheavy" begin

            for res in pdb
                heavy = findheavy(res)
                for i in 1:length(res.atoms)
                    if i in heavy
                        @test res.atoms[i].element != "H"
                    else
                        @test res.atoms[i].element == "H"
                    end
                end
            end
        end

        @testset "Interface between chain A and D" begin

            # Contacts: Chain A (4 residues) y Chain D (3 residues)
            C1 = @residuesdict pdb model "1" chain "A" group "ATOM" residue All
            C2 = @residuesdict pdb model "1" chain "D" group "ATOM" residue All

            TRUE = [ ("126",	"311"), ("183",	"309"), ("183",	"310"),
                     ("184",	"309"), ("187",	"309"), ("187",	"310"),
                     ("187",	"311"), ("187",	"312"), ("187",	"319"),
                     ("210",	"237"), ("211",	"312"), ("212",	"237"),
                     ("213",	"237"), ("213",	"312"), ("213",	"313"),
                     ("214",	"237") ]

            for (resnum1, resnum2) in TRUE
                res1 = C1[resnum1]
                res2 = C2[resnum2]

                @test contact(res1, res2, 6.5)

                if resnum2 == "309" && (resnum1 == "184" || resnum1 == "187")
                    @test  ionic(res1, res2)
                else
                    @test !ionic(res1, res2)
                end

                if (resnum1 == "211" && resnum2 == "312") || (resnum1 == "212" && resnum2 == "237")
                    @test  vanderwaals(res1, res2)
                else
                    @test !vanderwaals(res1, res2)
                end

                if resnum1 == "212" && resnum2 == "237"
                    @test  hydrophobic(res1, res2)
                else
                    @test !hydrophobic(res1, res2)
                end

                if resnum1 == "211" && resnum2 == "312"
                    @test  vanderwaalsclash(res1, res2)
                else
                    @test !vanderwaalsclash(res1, res2)
                end

                @test !aromaticsulphur(res1, res2)
                @test !pication(res1, res2)
                @test !disulphide(res1, res2)
                @test !aromatic(res1, res2)
                @test !hydrogenbond(res1, res2)
                @test !hydrogenbond(res1, res2)
                @test !covalent(res1, res2)
            end
        end

        @testset "Aromatic between chain A and B" begin

            C1 = @residuesdict pdb model "1" chain "A" group "ATOM" residue All
            C2 = @residuesdict pdb model "1" chain "B" group "ATOM" residue All

            @test aromatic(C1["36"],  C2["103"])
            @test aromatic(C1["94"],  C2["47"])
            @test aromatic(C1["94"],  C2["50"])
            @test aromatic(C1["96"],  C2["47"])
            @test aromatic(C1["98"],  C2["103"])
            @test aromatic(C1["135"], C2["174"])

            @test !aromatic(C1["96"],  C2["50"])
        end
    end

    @testset "1AKS" begin

        pdb = read(joinpath(pwd(), "data", "1AKS.xml.gz"), PDBML)

        CA = @residuesdict pdb model "1" chain "A" group "ATOM" residue All
        CB = @residuesdict pdb model "1" chain "B" group "ATOM" residue All

        @test  aromaticsulphur(CA["20"],  CB["157"])
        @test !aromaticsulphur(CA["29"],  CB["157"])

        @test  pication(CA["20"],  CB["159"])
        @test  pication(CA["40"],  CB["151"])
        @test  pication(CA["91"],  CB["234"])
        @test  pication(CA["91"],  CB["237"])
        @test !pication(CA["57"],  CB["215"])

        @test  disulphide(CA["22"],  CB["157"])
        @test  disulphide(CA["128"], CB["232"])
        @test  disulphide(CA["136"], CB["201"])

        @test  aromatic(CA["91"], CB["234"])
        @test  aromatic(CA["91"], CB["237"])
        @test !aromatic(CA["141"],CB["151"])
        @test !aromatic(CA["141"],CB["152"])
        @test !aromatic(CA["57"], CB["215"])
        @test !aromatic(CA["90"], CB["237"])

        @test  ionic(CA["87"], CB["245"]) # 1aksAB	A	87	LYS	NZ	B	245	ASN	OXT
        @test  ionic(CA["107"],CB["245"])
        @test  ionic(CA["135"],CB["159"])

        @test  covalent(CA["22"], CB["157"])
        @test  covalent(CA["128"],CB["232"])
        @test  covalent(CA["136"],CB["201"])

        for (res1, res2) in Tuple{String,String}[("136", "160"), ("17", "189"),
            ("20", "157"), ("140", "156"),
            ("138", "158"), ("144", "150"),
            ("137", "200")]

            @test hydrogenbond(CA[res1], CB[res2])
        end

        for (res1, res2) in Tuple{String,String}[("17", "221A"), ("138", "160"),
            ("87", "245"), ("89", "245"), ("16", "194"), ("17", "189"), ("17", "191"),
            ("17", "220"), ("20", "157"), ("22", "157"), ("124", "232"), ("128", "232"),
            ("134", "201"), ("136", "201"), ("137", "157"), ("143", "191"), ("16", "156"),
            ("124", "204"), ("124", "210"), ("129", "210"), ("143", "192"), ("144", "156"),
            ("47", "238"), ("47", "242"), ("48", "242"), ("51", "242"), ("53", "212"),
            ("53", "238"), ("55", "212"), ("103", "212"), ("103", "238"), ("105", "238"),
            ("105", "242"), ("123", "238"), ("16", "158"), ("21", "154"), ("22", "155"),
            ("27", "155"), ("30", "155"), ("47", "209"), ("53", "209"), ("72", "154"),
            ("121", "209"), ("123", "209"), ("138", "158"), ("141", "155"), ("20", "159"),
            ("135", "159"), ("137", "159"), ("99", "180"), ("100", "180"), ("29", "198"),
            ("141", "152"), ("144", "152"), ("51", "241"), ("89", "241"), ("100", "177"),
            ("103", "229"), ("105", "241"), ("89", "237"), ("91", "237"), ("99", "215"),
            ("103", "237"), ("105", "237"), ("101", "234"), ("103", "234"), ("138", "228"),
            ("143", "151"), ("29", "200"), ("121", "200"), ("123", "231"), ("123", "235"),
            ("124", "231"), ("124", "235"), ("129", "162"), ("134", "162"), ("136", "183"),
            ("136", "199"), ("138", "183"), ("138", "199"), ("138", "213"), ("51", "245")]

            @test hydrophobic(CA[res1], CB[res2])
        end

        for (res1, res2) in Tuple{String,String}[("136","160"), ("107","245"),
            ("16","194"), ("17","189"), ("124","232"), ("22","156"), ("142","192"),
            ("22","155"), ("140","155"), ("45","198"), ("57","195"), ("143","149"),
            ("140","194"), ("142","194"), ("19","157"), ("20","157"), ("22","157"),
            ("128","232"), ("134","201"), ("136","201"), ("138","157"), ("16","156"),
            ("124","210"), ("125","204"), ("127","210"), ("139","156"), ("140","156"),
            ("143","192"), ("55","196"), ("142","193"), ("103","212"), ("16","158"),
            ("53","209"), ("72","154"), ("122","209"), ("124","209"), ("138","158"),
            ("141","155"), ("136","159"), ("100","180"), ("29","198"), ("30","198"),
            ("134","161"), ("142","152"), ("144","152"), ("42","195"), ("43","195"),
            ("57","214"), ("72","153"), ("73","153"), ("74","153"), ("102","214"),
            ("143","150"), ("144","150"), ("145","146"), ("145","150"), ("145","147"), # OXT
            ("102","229"), ("91","237"), ("101","234"), ("143","151"), ("121","200"),
            ("134","162"), ("137","200"), ("100","177"), ("133","162")]

            @test vanderwaalsclash(CA[res1], CB[res2])
        end

        for (res1, res2) in Tuple{String,String}[("17","221A"), ("135","160"),
            ("89","245"), ("100","179"), ("101","179"), ("107","245"), ("134","202"),
            ("16","189"), ("16","194"), ("20","157"), ("22","157"), ("136","201"),
            ("137","157"), ("44","196"), ("40","193"), ("134","161"), ("124","209"),
            ("135","159"), ("44","198"), ("51","245"), ("102","214"), ("29","200"),
            ("17","189"), ("140","194"), ("141","194"), ("142","194"), ("19","157"),
            ("27","157"), ("124","232"), ("128","232"), ("134","201"), ("135","201"),
            ("138","157"), ("143","191"), ("16","156"), ("20","156"), ("21","156"),
            ("22","156"), ("122","204"), ("139","156"), ("140","156"), ("43","196"),
            ("124","204"), ("124","210"), ("125","204"), ("127","210"), ("129","210"),
            ("142","192"), ("143","192"), ("144","156"), ("17","188A"), ("18","188A"),
            ("44","197"), ("54","196"), ("55","196"), ("122","203"), ("142","193"),
            ("145","148"), ("47","238"), ("51","242"), ("137","158"), ("138","158"),
            ("53","212"), ("54","212"), ("55","212"), ("103","212"), ("103","238"),
            ("105","238"), ("123","238"), ("16","158"), ("136","159"), ("137","159"),
            ("21","154"), ("21","155"), ("22","155"), ("27","155"), ("30","155"),
            ("45","209"), ("47","209"), ("53","209"), ("45","198"), ("133","161"),
            ("71","154"), ("71","155"), ("72","154"), ("121","209"), ("122","209"),
            ("140","155"), ("141","154"), ("141","155"), ("17","188"), ("128","230"),
            ("98","180"), ("99","180"), ("100","180"), ("29","198"), ("30","198"),
            ("135","161"), ("138","198"), ("139","198"), ("141","152"), ("142","152"),
            ("144","152"), ("16","190"), ("17","146"), ("42","195"), ("72","153"),
            ("43","195"), ("57","195"), ("57","214"), ("58","195"), ("71","153"),
            ("132","164"), ("141","153"), ("143","149"), ("143","150"), ("144","150"),
            ("145","146"), ("145","147"), ("145","149"), ("73","153"), ("74","153"),
            ("145","150"), ("89","241"), ("100","177"), ("102","229"), ("105","241"),
            ("91","237"), ("92","237"), ("99","215"), ("138","228"), ("143","151"),
            ("103","237"), ("105","237"), ("91","234"), ("101","234"), ("103","234"),
            ("121","200"), ("124","231"), ("124","235"), ("132","162"), ("133","162"),
            ("134","162"), ("136","162"), ("136","199"), ("136","160"), ("138","160"),
            ("136","200"), ("137","199"), ("137","200"), ("138","183"), ("138","199")]

            @test vanderwaals(CA[res1], CB[res2])
        end
    end

    @testset "Vectorized contact/distance" begin

        code = "2VQC"
        pdb = read(txt(code), PDBFile)

        for criteria in ["All", "CA", "CB", "Heavy"]
            dist = distance(pdb, criteria=criteria)
            sq_d = squared_distance(pdb, criteria=criteria)
            cont = contact(pdb, 6.05, criteria=criteria)

            @test all(diag(cont))
            @test all(diag(dist) .== 0.0)
            @test all(diag(sq_d) .== 0.0)

            @test all((dist .<= 6.05) .== cont)
            @test all((sq_d .<= 36.6025) .== cont) # 6.05^2
        end
    end

    @testset "Proximity Mean" begin

        code = "2VQC"
        pdb = read(txt(code), PDBFile)
        residues = @residues pdb model "1" chain "A" group "ATOM" residue x -> x in ["62","64","65"]

        @test contact(residues, 6.05) == ( [1 1 0
                                            1 1 1
                                            0 1 1 ] .== 1 ) # All == Heavy (2VQC doesn't have Hs)

        @test proximitymean(residues,[1.0, 2.0, 3.0],6.05) == [2.0, 2.0, 2.0]
        @test proximitymean(residues,[10., 15., 30.],6.05) == [15., 20., 15.]

        @test proximitymean(residues,[1.0, 2.0, 3.0],6.05,include=true) == [3,12/3,5]./2.0
        @test proximitymean(residues,[10., 15., 30.],6.05,include=true) == [25,110/3,45]./2.0
    end
end
