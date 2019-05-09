@testset "SIFTS Mappings" begin

    sifts_file = joinpath(pwd(), "data", "2vqc.xml.gz")

    @testset "parse" begin

        dbs = [(dbUniProt, "P20220"), (dbPfam, "PF09645"), (dbNCBI, "244589"),
               (dbPDB, "2vqc"), (dbSCOP, "153426"), (dbPDBe, "2vqc")]
        for (to,toid) in dbs, (from, fromid) in dbs
            map = siftsmapping(sifts_file, from, fromid, to, toid)
            for (k,v) in map
                to == from && @test k == v
            end
        end
    end

    @testset "missings = false" begin

        map = siftsmapping(sifts_file,dbPDBe,"2vqc",dbPDB,"2vqc",chain="A",missings=false)
        @test_throws KeyError map["9"]  # Missing
        @test_throws KeyError map["80"] # Missing
        @test_throws KeyError map["1"]  # Missing
        @test map["10"] == "4"
        @test map["79"] == "73"
    end

    @testset "missings = true" begin

        map = siftsmapping(sifts_file, dbPDBe, "2vqc", dbPDB, "2vqc", chain="A")
        @test map["9"]  == "3"   # Missing
        @test map["80"] == "74"  # Missing
        @test map["1"]  == "-5"  # Missing # Negative Resnum
        @test map["10"] == "4"
        @test map["79"] == "73"
    end

    @testset "Insert codes" begin
    # 1SSX : Residues with insert codes: 15A 15B

        _1ssx_file = joinpath(pwd(), "data", "1ssx.xml.gz")

        map = siftsmapping(_1ssx_file, dbPDBe, "1ssx", dbPDB, "1ssx", chain="A")
        residue_A = map["1"]
        residue_B = map["2"]
        residue_C = map["3"]
        @test residue_A == "15A"
        @test residue_B == "15B"
        @test residue_C == "16"

        mapII = siftsmapping(_1ssx_file, dbPDB, "1ssx", dbUniProt, "P00778", chain="A")
        @test mapII["15A"] == "200"
        @test mapII["15B"] == "201"
        @test mapII["16"]  == "202"
    end

    @testset "Multiple InterProt annotations" begin
    # 1CBN : Multiple InterProt annotations, the last is used.
    # Identical PDBe ResNum for Residue 22:
    #
    #         <residue dbSource="PDBe" dbCoordSys="PDBe" dbResNum="22" dbResName="SER">
    #           <crossRefDb dbSource="PDB" dbCoordSys="PDBresnum" dbAccessionId="1cbn" dbResNum="22" dbResName="SER" dbChainId="A"/>
    #           <crossRefDb dbSource="UniProt" dbCoordSys="UniProt" dbAccessionId="P01542" dbResNum="22" dbResName="P"/>
    #           ....
    #         <residue dbSource="PDBe" dbCoordSys="PDBe" dbResNum="22" dbResName="PRO">
    #           <crossRefDb dbSource="PDB" dbCoordSys="PDBresnum" dbAccessionId="1cbn" dbResNum="22" dbResName="PRO" dbChainId="A"/>
    #           <crossRefDb dbSource="UniProt" dbCoordSys="UniProt" dbAccessionId="P01542" dbResNum="22" dbResName="P"/>
    #           ...

        _1cbn_file = joinpath(pwd(), "data", "1cbn.xml.gz")

        map = siftsmapping(_1cbn_file, dbPDBe, "1cbn", dbInterPro, "IPR001010", chain="A")
        @test_throws KeyError map["1"] # Without InterPro
        @test map["2"] == "2" # Same ResNum for different InterPros

        mapII = read(_1cbn_file, SIFTSXML, chain="A")
        @test length(mapII[1].InterPro) == 0 # Without InterPro
        @test length(mapII[2].InterPro) == 4 # Same ResNum for different InterPros
    end

    @testset "None" begin
    # 4CPA : It has "None"
    # <residue dbSource="PDBe" dbCoordSys="PDBe" dbResNum="2" dbResName="GLX">
    #   <crossRefDb dbSource="PDB" dbCoordSys="PDBresnum" dbAccessionId="4cpa" dbResNum="2" dbResName="GLX" dbChainId="J"/>
    #   <crossRefDb dbSource="SCOP" dbCoordSys="PDBresnum" dbAccessionId="118683" dbResNum="2" dbResName="GLX" dbChainId="J"/>
    #   <crossRefDb dbSource="InterPro" dbCoordSys="UniProt" dbAccessionId="IPR011052" dbResNum="None" dbResName="None" dbEvidence="SSF57027"/>
    #   <crossRefDb dbSource="InterPro" dbCoordSys="UniProt" dbAccessionId="IPR021142" dbResNum="None" dbResName="None" dbEvidence="PD884054"/>
    #   <crossRefDb dbSource="InterPro" dbCoordSys="UniProt" dbAccessionId="IPR004231" dbResNum="None" dbResName="None" dbEvidence="PF02977"/>
    #   <residueDetail dbSource="PDBe" property="codeSecondaryStructure">T</residueDetail>
    #   <residueDetail dbSource="PDBe" property="nameSecondaryStructure">loop</residueDetail>
    # </residue>

        _4cpa_file = joinpath(pwd(), "data", "4cpa.xml.gz")

        map = read(_4cpa_file, SIFTSXML, chain="J")
        res = filter(r -> r.PDBe.number == "2", map)[1]
        @test length(res.InterPro) == 3
        for i in 1:3
            @test res.InterPro[i].name == ""
            @test res.InterPro[i].number == ""
        end
    end

    @testset "NMR" begin
    # 1AS5 : NMR

        _1as5_file = joinpath(pwd(), "data", "1as5.xml.gz")

        map = siftsmapping(_1as5_file, dbPDBe,"1as5", dbUniProt, "P56529", chain="A")
                          # missings=true : NMR there are not missing residues
        @test map["23"] == "73"
        @test_throws KeyError map["24"] # Without UniProt

        mapII = siftsmapping(_1as5_file, dbPDBe,"1as5", dbPDB, "1as5", chain="A")
        @test mapII["24"] == "24"
    end

    @testset "Inserted residues lack insertion code" begin
    # 1DPO : Inserted residues lack insertion letters
    # Single unnamed chain in 1DPO contains insertions at postions 184 (Gly, Phe),
    # 188 (Gly, Lys), and 221 (Ala, Leu) but no insertion letters.

        _1dpo_file = joinpath(pwd(), "data", "1dpo.xml.gz")
        map = siftsmapping(_1dpo_file, dbPDBe, "1dpo", dbPDB, "1dpo", chain="A")
        # Unnamed chain is "A" in SIFTS
        @test map["164"] == "184"
        @test map["165"] == "184A" # Has insertion code in SIFTS
        @test map["169"] == "188"
        @test map["170"] == "188A" # Has insertion code in SIFTS
        @test map["198"] == "221"
        @test map["199"] == "221A" # Has insertion code in SIFTS
    end

    @testset "Insertion block" begin
    # 1IGY : Insertions have more than one copy of the same amino acid in a single
    # insertion block. For example, chain B in 1IGY contains a block of four residues
    # inserted at sequence position 82. The block contains Leu-Ser-Ser-Leu.

        _1igy_file = joinpath(pwd(), "data", "1igy.xml.gz")
        map = siftsmapping(_1igy_file, dbPDBe, "1igy", dbCATH, "2.60.40.10", chain="B")
        @test map["82"] == "82"
        @test map["83"] == "82A"
        @test map["84"] == "82B"
        @test map["85"] == "82C"
    end

    @testset "1HAG" begin
    # 1HAG : Chain E begins with 1H, 1G, 1F, ... 1A, then 1 (in reverse alphabetic order)

        _1hag_file = joinpath(pwd(), "data", "1hag.xml.gz")
        map = siftsmapping(_1hag_file, dbPDBe, "1hag", dbPDB, "1hag", chain="E")
        @test map["1"] == "1H"
        @test map["2"] == "1G"
        @test map["3"] == "1F"
        @test map["4"] == "1E"
        @test map["5"] == "1D"
        @test map["6"] == "1C"
        @test map["7"] == "1B"
        @test map["8"] == "1A"
        @test map["9"] == "1"
    end
end

@testset "1NSA" begin
# 1NSA : Contains a single (unnamed) protein chain with sequence 7A-95A that continues 4-308.

    _1nsa_file = joinpath(pwd(), "data", "1nsa.xml.gz")
    mapping = read(_1nsa_file, SIFTSXML)

    @testset "findall & read" begin

        four = findall(db -> db.id == "1nsa" && db.number == "4" , mapping, dbPDB)[1]
        @test findall(db -> db.id == "1nsa" && db.number == "95A", mapping, dbPDB)[1]+1 == four
        @test mapping[findall(db -> db.id == "1nsa" && db.number == "95A",mapping,dbPDB)][1].PDB.number == "95A"
    end

    @testset "capture fields" begin

        cap = map(mapping) do res
            number = get(res, dbPDB, :number, "")
            if get(res, dbPDB, :id, "") == "1nsa" && number == "95A"
                number
            else
                missing
            end
        end
        @test cap[ [!ismissing(x) for x in cap] ][1] == "95A"
    end
end

@testset "find & filter" begin
# 1IAO : Contains in chain B (in this order) 1S, 323P-334P, 6-94, 94A, 95-188, 1T, 2T

    mapp = read(joinpath(pwd(), "data", "1iao.xml.gz"), SIFTSXML)

    @test filter(db -> db.id == "1iao" && db.number == "1S" && db.chain == "B", mapp, dbPDB)[1].PDBe.number == "1"
    i = findall(db -> db.id == "1iao" && db.number == "1S" && db.chain == "B", mapp, dbPDB)[1]
    res = mapp[i+2]
    @test get(res,dbPDB,:id,"") == "1iao" &&  get(res,dbPDB,:chain,"") == "B" &&
          get(res,dbPDB,:number,"") == "323P"
end

@testset "MIToS 1.0 error" begin

    sf = joinpath(pwd(), "data", "4gcr.xml.gz")

    mapping = siftsmapping(sf, dbPfam, "PF00030", dbPDB, "4gcr")
    @test mapping["3"] == "2"
end

@testset "download" begin

    pdb = "2vqc"
    mapping = read(joinpath(pwd(), "data", "$(pdb).xml.gz"), SIFTSXML)

    @test_throws AssertionError downloadsifts(pdb, filename="bad_name.txt")
    @test_throws ErrorException downloadsifts("2vqc_A")
    # FTP causes problem with Travis-CI
    # filename = downloadsifts(pdb, filename=tempname()*".xml.gz")
    # @test length(read(filename, SIFTSXML)) == length(mapping)
end

@testset "Ensembl" begin
# 18GS has multiple Ensembl annotations for each residue
# 18GS also has EC and GO annotations

    mapping = read(joinpath(pwd(), "data", "18gs.xml.gz"), SIFTSXML)

    last_res = mapping[end]
    @test length(last_res.Ensembl) == 2
    @test last_res.Ensembl[1].id == last_res.Ensembl[2].id
    @test last_res.Ensembl[2].transcript == "ENST00000642444"
    @test last_res.Ensembl[1].transcript == "ENST00000398606"
    @test last_res.Ensembl[2].transcript == "ENST00000642444"
    @test last_res.Ensembl[1].translation == "ENSP00000381607"
    @test last_res.Ensembl[2].translation == "ENSP00000493538"
    @test last_res.Ensembl[1].exon == "ENSE00001921020"
    @test last_res.Ensembl[2].exon == "ENSE00003822108"
end
