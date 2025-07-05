let
    sifts_file = joinpath(@__DIR__, "..", "..", "test", "data", "18gs.xml.gz")
    xdoc = Utils._get_xml_document(sifts_file)
    residue = first(
        SIFTS._get_residues(first(SIFTS._get_segments(first(SIFTS._get_entities(xdoc))))),
    )
    missing_residue, sscode, ssname = SIFTS._get_details(residue)
    SUITE["SIFTS"]["SIFTSResidue"]["18gs"] =
        @benchmarkable SIFTS.SIFTSResidue($residue, $missing_residue, $sscode, $ssname)
    SIFTS.LightXML.free(xdoc)
end
