let
    cif_file = joinpath(@__DIR__, "..", "..", "test", "data", "2vqc.cif")
    SUITE["PDB"]["read"]["MMCIFFile"] = @benchmarkable read_file($cif_file, MMCIFFile)
    residues = read_file(cif_file, MMCIFFile)
    SUITE["PDB"]["_pdbresidues_to_mmcifdict"]["2vqc"] =
        @benchmarkable PDB._pdbresidues_to_mmcifdict($residues)
end
