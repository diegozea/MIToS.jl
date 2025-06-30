let
    pdb_file = joinpath(@__DIR__, "..", "..", "test", "data", "2hhb.pdb.gz")
    hemoglobin = read_file(pdb_file, PDBFile, group = "ATOM", model = "1")
    chain_a = select_residues(hemoglobin, model = "1", chain = "A", group = "ATOM")
    chain_c = select_residues(hemoglobin, model = "1", chain = "C", group = "ATOM")
    matches = zip(eachindex(chain_a), eachindex(chain_c))
    SUITE["PDB"]["_get_matched_Cαs"]["hemoglobin"] =
        @benchmarkable PDB._get_matched_Cαs($chain_a, $chain_c, $matches)
end
