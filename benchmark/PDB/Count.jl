let
    pdb_file = joinpath(@__DIR__, "..", "..", "test", "data", "1CBN.pdb")
    struc = read_file(pdb_file, PDBFile)

    SUITE["PDB"]["count_alanine"]["1CBN"] =
        @benchmarkable count(res -> res.id.name == "ALA", $struc)
end
