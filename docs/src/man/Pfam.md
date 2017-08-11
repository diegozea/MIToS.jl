
# Pfam

MIToS defines methods and types useful for any MSA. The `Pfam` module uses other MIToS modules in the context of Pfam MSAs, where it’s possible to us determine how structure and sequence information should be mapped. This module defines functions that go from a Pfam MSA to the protein contact prediction performance of pairwise scores estimated from that MSA.

## Features

- [**Download and read**](#Getting-a-Pfam-MSA) Pfam MSAs
- Obtain [**PDB information**](#Getting-PDB-information-from-an-MSA) from alignment annotations
- [**Map**](#Getting-PDB-information-from-an-MSA) between sequence/alignment residues/columns and PDB structures
- Measure of [**AUC**](#PDB-contacts-and-AUC) (ROC curve) for [**protein contact**](#PDB-contacts-and-AUC) prediction of MI scores


```julia
using MIToS.Pfam
```

## Contents

- [Getting a Pfam MSA](#Getting-a-Pfam-MSA)
- [Getting PDB information from and MSA](#Getting-PDB-information-from-an-MSA)
- [PDB contacts and AUC](#PDB-contacts-and-AUC)


```julia
# Truncate IJulia outputs at:
ENV["LINES"]   = 15 
ENV["COLUMNS"] = 60;
```

<a href="#"><i class="fa fa-arrow-up"></i></a>

## Getting a Pfam MSA

The function `downloadpfam` takes a Pfam accession and downloads a Pfam MSA in Stockholm format. Use `read` function and the `Stockholm` `Format` to get a `AnnotatedMultipleSequenceAlignment` object with the MSA and its Pfam annotations. You must set `generatemapping` and `useidcoordinates` to `true` the first time you read the downloaded MSA. This is necessary to some of the methods in the `Pfam` module.


```julia
pfamfile = downloadpfam("PF12464")
```

      % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                     Dload  Upload   Total   Spent    Left  Speed
    100 73781  100 73781    0     0  32697      0  0:00:02  0:00:02 --:--:-- 32704





    "PF12464.stockholm.gz"




```julia
msa = read(pfamfile, Stockholm, generatemapping=true, useidcoordinates=true)
```




    1157x53 MIToS.MSA.AnnotatedMultipleSequenceAlignment:
     R  M  S  S  G  Q  L  Y  I  …  M  K  G  F  F  G  A  C  G
     -  C  R  L  G  E  L  Y  N     I  R  D  L  L  G  K  T  G
     R  M  L  A  G  L  P  Y  R     L  A  E  I  L  G  K  C  G
     K  M  L  A  G  E  L  Y  D     L  K  E  L  L  G  -  -  -
     R  M  L  L  G  L  P  Y  K     T  R  E  I  L  Y  K  V  G
     ⋮              ⋮           ⋱     ⋮              ⋮      
     K  M  T  A  G  E  W  Y  C     L  A  L  L  F  A  K  -  -
     K  M  I  A  G  D  L  Y  F     V  K  E  T  F  G  S  V  G
     -  -  -  -  -  L  P  Y  Y     L  R  K  L  L  G  K  T  G
     R  M  I  S  G  M  L  Y  N  …  L  R  E  I  L  G  S  I  -
     K  M  I  Q  G  E  L  Y  Y     -  -  -  -  -  -  -  -  -



<a href="#"><i class="fa fa-arrow-up"></i></a>

## Getting PDB information from an MSA

The function `getseq2pdb` parses the MSA annotations to return a `Dict` from the sequence identifier in the MSA to PDB and chain codes.


```julia
getseq2pdb(msa)
```




    Dict{ASCIIString,Array{Tuple{ASCIIString,ASCIIString},1}} with 5 entries:
      "Q18A66_PEPD6/6-5… => [("4ISX","A"),("3SRT","B"),("3SRT",…
      "THGA_ECOLI/8-59"  => [("1KRR","B"),("1KRU","A"),("1KRV",…
      "MAA_ECOLI/7-58"   => [("1OCX","C"),("1OCX","A"),("1OCX",…
      "Q9KLB0_VIBCH/8-5… => [("3NZ2","J"),("3NZ2","D"),("3NZ2",…
      "Q81N16_BACAN/7-5… => [("3HJJ","C"),("3IGJ","C"),("3IGJ",…



Once you know the association between PDB chains and sequences, you can use that information together with the `msacolumn2pdbresidue` function to get the PDB residue number that correspond to each MSA column for given a determined sequence and PDB chain. That function downloads information from SIFTS to generate the mapping.


```julia
col2res = msacolumn2pdbresidue(msa, "MAA_ECOLI/7-58", "1OCX", "C")
```

      % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                     Dload  Upload   Total   Spent    Left  Speed
    100 22683  100 22683    0     0   4074      0  0:00:05  0:00:05 --:--:--  4436





    Dict{Int64,ASCIIString} with 51 entries:
      46  => "21"
      134 => "48"
      136 => "50"
      55  => "29"
      66  => "38"
      58  => "32"
      59  => "33"
      142 => "56"
      139 => "53"
      57  => "31"
      ⋮   => ⋮



The returned dictionary can be used to get the PDB residue associated to each column (using the `msaresidues` function)...


```julia
using MIToS.PDB
pdbfile = downloadpdb("1OCX")
pdb = read(pdbfile, PDBML)
resdict = @residuesdict pdb model "1" chain "C" group "ATOM" residue "*"

msaresidues(msa, resdict, col2res)
```

      % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                     Dload  Upload   Total   Spent    Left  Speed
      0     0    0     0    0     0      0      0 --:--:-- --:--:-- --:--:--     0
    100  355k  100  355k    0     0   178k      0  0:00:01  0:00:01 --:--:--  757k





    DataStructures.OrderedDict{Int64,MIToS.PDB.PDBResidue} with 51 entries:
      14 => PDBResidue:…
      15 => PDBResidue:…
      16 => PDBResidue:…
      17 => PDBResidue:…
      18 => PDBResidue:…
      19 => PDBResidue:…
      20 => PDBResidue:…
      21 => PDBResidue:…
      22 => PDBResidue:…
      24 => PDBResidue:…
      ⋮  => ⋮



...or to delete the columns without PDB residues (using the `hasresidues` function):


```julia
using MIToS.MSA
filtercolumns!(msa, hasresidues(msa, col2res))
```




    1157x51 MIToS.MSA.AnnotatedMultipleSequenceAlignment:
     R  M  S  S  G  Q  L  Y  I  …  L  M  K  G  F  F  G  A  C
     -  C  R  L  G  E  L  Y  N     I  I  R  D  L  L  G  K  T
     R  M  L  A  G  L  P  Y  R     L  L  A  E  I  L  G  K  C
     K  M  L  A  G  E  L  Y  D     I  L  K  E  L  L  G  -  -
     R  M  L  L  G  L  P  Y  K     L  T  R  E  I  L  Y  K  V
     ⋮              ⋮           ⋱           ⋮              ⋮
     K  M  T  A  G  E  W  Y  C     A  L  A  L  L  F  A  K  -
     K  M  I  A  G  D  L  Y  F     R  V  K  E  T  F  G  S  V
     -  -  -  -  -  L  P  Y  Y     R  L  R  K  L  L  G  K  T
     R  M  I  S  G  M  L  Y  N  …  M  L  R  E  I  L  G  S  I
     K  M  I  Q  G  E  L  Y  Y     -  -  -  -  -  -  -  -  -



<a href="#"><i class="fa fa-arrow-up"></i></a>

### PDB contacts and AUC

The `Dict` between MSA columns and PDB residue number also can be used to generate a protein contact map associated to the MSA.


```julia
cmap = msacontacts(msa, resdict, col2res)
```




    51x51 PairwiseListMatrices.PairwiseListMatrix{Float64,false}:
     NaN      1.0    1.0    1.0  …    0.0    0.0    0.0
       1.0  NaN      1.0    1.0       0.0    0.0    0.0
       1.0    1.0  NaN      1.0       0.0    0.0    0.0
       1.0    1.0    1.0  NaN         0.0    0.0    0.0
       1.0    1.0    1.0    1.0       0.0    0.0    0.0
       ⋮                         ⋱                  ⋮  
       0.0    0.0    0.0    0.0       1.0    0.0    0.0
       0.0    0.0    0.0    0.0       1.0    1.0    1.0
       0.0    0.0    0.0    0.0     NaN      1.0    1.0
       0.0    0.0    0.0    0.0       1.0  NaN      1.0
       0.0    0.0    0.0    0.0  …    1.0    1.0  NaN  



That protein contact map can be used to calculate the Area Under the ROC Curve for a given score with the `AUC` function.


```julia
using MIToS.Information
ZMIp, MIp = buslje09(msa)

AUC(ZMIp, cmap)
```




    0.7549031963326455


