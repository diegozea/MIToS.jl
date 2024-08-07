"""
The `Pfam` module, defines functions to measure the protein contact prediction performance of information measure between column pairs from a Pfam MSA.

**Features**

  - Read and download Pfam MSAs
  - Obtain PDB information from alignment annotations
  - Map between sequence/alignment residues/columns and PDB structures
  - Measure of AUC (ROC curve) for contact prediction of MI scores

```julia
using MIToS.Pfam
```
"""
module Pfam

using MIToS.Utils
using MIToS.MSA
using MIToS.SIFTS
using MIToS.PDB
using MIToS.Information
using PairwiseListMatrices
using NamedArrays
using OrderedCollections

export  # Download
    downloadpfam,
    # PDB
    Stockholm,
    getseq2pdb,
    msacolumn2pdbresidue,
    hasresidues,
    msacontacts,
    msaresidues,
    getcontactmasks,
    # Utils
    read_file,
    parse_file,
    write_file,
    print_file

include("Download.jl")
include("PDB.jl")

end
