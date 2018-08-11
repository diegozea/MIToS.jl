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
using DataStructures
using Requires

export  # Download
        downloadpfam,
        # PDB
        Stockholm,
        getseq2pdb,
        msacolumn2pdbresidue,
        hasresidues,
        msacontacts,
        msaresidues,
        getcontactmasks

include("Download.jl")
include("PDB.jl")

@require ROCAnalysis begin

    """
    `AUC(scores_list::Vector, true_contacts::BitVector, false_contacts::BitVector)`

    Returns the Area Under a ROC (Receiver Operating Characteristic) Curve (AUC) of the
    `scores_list` for `true_contacts` prediction. The three vectors should have the same
    length and `false_contacts` should be `true` where there are not contacts. You need to
    do `using ROCAnalysis` before using this function.
    """
    function ROCAnalysis.AUC(scores_list::Vector{T},
                             true_contacts::BitVector,
                             false_contacts::BitVector) where T
        ROCAnalysis.AUC(ROCAnalysis.roc(
            scores_list[true_contacts  .& .!isnan.(scores_list)],
            scores_list[false_contacts .& .!isnan.(scores_list)]))
    end

    """
    `AUC(scores::PairwiseListMatrix, true_contacts::BitVector, false_contacts::BitVector)`

    Returns the Area Under a ROC (Receiver Operating Characteristic) Curve (AUC) of the
    `scores` for `true_contacts` prediction. `scores`, `true_contacts` and `false_contacts`
    should have the same number of elements and `false_contacts` should be `true` where
    there are not contacts. You need to do `using ROCAnalysis` before using this function.
    """
    function ROCAnalysis.AUC(scores::NamedArray{L,2,PairwiseListMatrix{L,false,VL},NL},
                             true_contacts::BitVector,
                             false_contacts::BitVector) where {L,VL,NL}
        ROCAnalysis.AUC(getlist(getarray(scores)), true_contacts, false_contacts)
    end

    """
    `AUC(scores::PairwiseListMatrix, msacontacts::PairwiseListMatrix)`

    Returns the Area Under a ROC (Receiver Operating Characteristic) Curve (AUC) of the
    `scores` for `msacontact` prediction. `score` and `msacontact` lists are vinculated
    (inner join) by their labels (i.e. column number in the file). `msacontact` should have
    1.0 for true contacts and 0.0 for not contacts (NaN or other numbers for missing
    values). You need to do `using ROCAnalysis` before using this function.
    """
    function ROCAnalysis.AUC(
    scores::NamedArray{L,2,PairwiseListMatrix{L,false,VL},NL},
    msacontacts::NamedArray{R,2,PairwiseListMatrix{L,false,VR},NR}) where {L <: AbstractFloat, R <: AbstractFloat,VL,VR,NL,NR}
        sco, con = join(scores, msacontacts, kind=:inner)
        true_contacts, false_contacts = getcontactmasks(con)
        ROCAnalysis.AUC(sco, true_contacts, false_contacts)
    end

    export AUC

end

end
