let
    pdb_file = joinpath(@__DIR__, "..", "..", "test", "data", "1CBN.pdb")
    struc = read_file(pdb_file, PDBFile, model = "1", chain = "A", group = "ATOM")

    SUITE["PDB"]["distance"]["1CBN_20_30"] = @benchmarkable distance($struc[20], $struc[30])
end
