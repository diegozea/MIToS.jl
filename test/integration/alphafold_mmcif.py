"""Check MIToS exports with an installed, unmodified AlphaFold mmCIF parser.

Run via AlphaFoldMMCIF.jl. No AlphaFold inference, weights, template search, or
ColabFold installation is needed. AlphaFold must be importable on PYTHONPATH.
"""

import sys
from pathlib import Path

import Bio
from alphafold.data import mmcif_parsing


def parse(directory, name):
    result = mmcif_parsing.parse(
        file_id=name,
        mmcif_string=(directory / f"{name}.cif").read_text(),
        catch_all_errors=False,
    )
    assert not result.errors, result.errors
    assert result.mmcif_object is not None
    obj = result.mmcif_object
    for chain, sequence in obj.chain_to_seqres.items():
        mapping = obj.seqres_to_structure[chain]
        assert sorted(mapping) == list(range(len(sequence)))
        for index, residue in mapping.items():
            assert not residue.is_missing
            position = residue.position
            # Verify the sequence-to-coordinate map actually retrieves each residue.
            atom_residue = obj.structure[chain][
                residue.hetflag, position.residue_number, position.insertion_code
            ]
            assert atom_residue.resname == residue.name
    assert "release_date" not in obj.header  # No fabricated date.
    assert obj.header["structure_method"] == "?"
    return obj


directory = Path(sys.argv[1])
inserted = parse(directory, "insertions")
renumbered = parse(directory, "renumbered")
assert inserted.chain_to_seqres == renumbered.chain_to_seqres
positions = inserted.seqres_to_structure["A"]
assert positions[0].position.residue_number == 15
assert positions[0].position.insertion_code == "A"
assert positions[1].position.residue_number == 15
assert positions[1].position.insertion_code == "B"
for index, residue in renumbered.seqres_to_structure["A"].items():
    assert residue.position.residue_number == index + 1
    assert residue.position.insertion_code == " "
nmr = parse(directory, "nmr")
assert len(nmr.chain_to_seqres["A"]) == 24
assert set(nmr.raw_string["_atom_site.pdbx_PDB_model_num"]) == {"1"}
modified = parse(directory, "modified")
assert modified.chain_to_seqres == {"AA": "AMG", "Z": "AMG"}
assert modified.seqres_to_structure["AA"][1].hetflag == "H_MSE"
print(f"AlphaFold parser integration passed for 4 exports (Biopython {Bio.__version__})")
