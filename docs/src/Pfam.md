```@setup log
@info "Pfam docs"
```

# [Pfam](@id Module-Pfam)

MIToS defines methods and types useful for any MSA. The `Pfam` module uses other MIToS
modules in the context of Pfam MSAs, where it’s possible to us determine how structure and
sequence information should be mapped. This module defines functions that go from a Pfam
MSA to the protein contact prediction performance of pairwise scores estimated from that MSA.

```julia
using MIToS.Pfam # to load the Pfam module
```  

## Features

- [**Download and read**](@ref Getting-a-Pfam-MSA) Pfam MSAs.
- Obtain [**PDB information**](@ref Getting-PDB-information-from-an-MSA) from alignment annotations.
- [**Map**](@ref Getting-PDB-information-from-an-MSA) between sequence/alignment residues/columns and PDB structures.
- Measure of [**AUC**](@ref PDB-contacts-and-AUC) (ROC curve) for [**protein contact**](@ref PDB-contacts-and-AUC) prediction of MI scores.

## Contents

```@contents
Pages = ["Pfam.md"]
Depth = 4
```

## [Getting a Pfam MSA](@id Getting-a-Pfam-MSA)

The function `downloadpfam` takes a Pfam accession and downloads a Pfam MSA in Stockholm
format. Use `read` function and the `Stockholm` `FileFormat` to get a
`AnnotatedMultipleSequenceAlignment` object with the MSA and its Pfam annotations.
You must set `generatemapping` and `useidcoordinates` to `true` the first time you read
the downloaded MSA. This is necessary to some of the methods in the `Pfam` module.  

```@example pfam_example
using MIToS.Pfam
pfamfile = downloadpfam("PF12464")
msa = read(pfamfile, Stockholm, generatemapping=true, useidcoordinates=true)
```

## [Getting PDB information from an MSA](@id Getting-PDB-information-from-an-MSA)  

The function `getseq2pdb` parses the MSA annotations to return a `Dict` from the sequence
identifier in the MSA to PDB and chain codes.  

```@example pfam_example
getseq2pdb(msa)
```

Once you know the association between PDB chains and sequences, you can use that
information together with the `msacolumn2pdbresidue` function to get the PDB residue
number that correspond to each MSA column for given a determined sequence and PDB chain.
That function downloads information from SIFTS to generate the mapping.  

```@example pfam_example
col2res = msacolumn2pdbresidue(msa, "MAA_ECOLI/7-58", "1OCX", "C")
```

The returned dictionary can be used to get the PDB residue associated to each column
(using the `msaresidues` function)...  

```@example pfam_example
using MIToS.PDB
pdbfile = downloadpdb("1OCX")
pdb = read(pdbfile, PDBML)
resdict = @residuesdict pdb model "1" chain "C" group "ATOM" residue All

msaresidues(msa, resdict, col2res)
```

...or to delete the columns without PDB residues (using the `hasresidues` function):  

```@example pfam_example
using MIToS.MSA
filtercolumns!(msa, hasresidues(msa, col2res))
```

### [PDB contacts and AUC](@id PDB-contacts-and-AUC)  

The `Dict` between MSA columns and PDB residue number also can be used to generate a
protein contact map associated to the MSA.  

```@example pfam_example
cmap = msacontacts(msa, resdict, col2res)
```

That protein contact map can be used to calculate the Area Under the ROC Curve for a given
score with the `AUC` function.  

```@example pfam_example
using MIToS.Information
ZMIp, MIp = buslje09(msa)

using ROCAnalysis # You need to load ROCAnalysis to use the AUC function

AUC(ZMIp, cmap)
```
