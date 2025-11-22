let
    sifts_file = joinpath(@__DIR__, "..", "..", "test", "data", "18gs.xml.gz")

    global function build_residue()
        xdoc = Utils._get_xml_document(sifts_file)
        residue = first(
            SIFTS._get_residues(
                first(SIFTS._get_segments(first(SIFTS._get_entities(xdoc)))),
            ),
        )
        return residue, xdoc
    end

    residue_for_constructor, xdoc_for_constructor = build_residue()
    missing_residue, sscode, ssname = SIFTS._get_details(residue_for_constructor)
    SUITE["SIFTS"]["SIFTSResidue"]["18gs"] = @benchmarkable SIFTS.SIFTSResidue(
        $residue_for_constructor,
        $missing_residue,
        $sscode,
        $ssname,
    )
    SUITE["SIFTS"]["ResidueDetails"]["_get_details"] =
        @benchmarkable SIFTS._get_details(residue) setup=((residue, xdoc) = build_residue()) teardown=(SIFTS.LightXML.free(
            xdoc,
        ))
    SUITE["SIFTS"]["ResidueDetails"]["_is_missing"] =
        @benchmarkable SIFTS._is_missing(residue) setup=((residue, xdoc) = build_residue()) teardown=(SIFTS.LightXML.free(
            xdoc,
        ))
    SIFTS.LightXML.free(xdoc_for_constructor)
end
