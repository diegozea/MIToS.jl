# Optional integration check; not part of Pkg.test and adds no Julia dependencies.
# Put an AlphaFold checkout on PYTHONPATH and install its mmcif_parsing dependencies.
# ALPHAFOLD_PYTHON=/path/to/python julia --project test/integration/AlphaFoldMMCIF.jl
using MIToS.PDB
import BioStructures

const DATA = joinpath(@__DIR__, "..", "data")

mktempdir() do directory
    inserted = read_file(joinpath(DATA, "1SSX.pdb"), PDBFile; chain = "A", group = "ATOM")
    for (name, renumber) in (("insertions", false), ("renumbered", true))
        cif = alphafold_mmcifdict(inserted; entry_id = name, renumber_auth = renumber)
        BioStructures.writemmcif(joinpath(directory, "$name.cif"), cif)
    end

    nmr = read_file(joinpath(DATA, "1AS5.pdb"), PDBFile)
    filter!(r -> r.id.name != "NH2", nmr) # Exclude the non-amino-acid terminal cap.
    cif = alphafold_mmcifdict(
        nmr;
        entry_id = "nmr",
        model = "14",
        chem_comp_types = Dict("HYP" => "L-peptide linking"),
    )
    BioStructures.writemmcif(joinpath(directory, "nmr.cif"), cif)

    # Repeat a sequence in a multi-character chain, with a modified HETATM monomer.
    residues = PDBResidue[]
    for chain in ("AA", "Z"), (index, name) in enumerate(("ALA", "MSE", "GLY"))
        push!(
            residues,
            PDBResidue(
                PDBResidueIdentifier(
                    "",
                    string(index + 20),
                    name,
                    name == "MSE" ? "HETATM" : "ATOM",
                    "5",
                    chain,
                ),
                [PDBAtom(Coordinates(index, 0, 0), "CA", "C", 1.0, "0.0", "", "")],
            ),
        )
    end
    cif = alphafold_mmcifdict(
        residues;
        entry_id = "modified",
        chem_comp_types = Dict("MSE" => "L-peptide linking"),
    )
    BioStructures.writemmcif(joinpath(directory, "modified.cif"), cif)

    python = get(ENV, "ALPHAFOLD_PYTHON", "python")
    checker = joinpath(@__DIR__, "alphafold_mmcif.py")
    run(`$python $checker $directory`)
end
