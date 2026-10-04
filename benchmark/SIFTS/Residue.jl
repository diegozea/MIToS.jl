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

    SUITE["SIFTS"]["SIFTSResidue"]["18gs"] =
        @benchmarkable SIFTS.SIFTSResidue(residue, missing_residue, sscode, ssname) setup =
            (
                (residue, xdoc) = build_residue();
                (missing_residue, sscode, ssname) = SIFTS._get_details(residue)
            ) teardown = (SIFTS.LightXML.free(xdoc))
    SUITE["SIFTS"]["ResidueDetails"]["_get_details"] =
        @benchmarkable SIFTS._get_details(residue) setup=((residue, xdoc) = build_residue()) teardown=(SIFTS.LightXML.free(
            xdoc,
        ))
    SUITE["SIFTS"]["ResidueDetails"]["_is_missing"] =
        @benchmarkable SIFTS._is_missing(residue) setup=((residue, xdoc) = build_residue()) teardown=(SIFTS.LightXML.free(
            xdoc,
        ))
end
