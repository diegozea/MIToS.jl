
# SIFTS

The `SIFTS` module of MIToS allows to obtain the residue-level mapping between databases
stored in the SIFTS XML files. It makes easy to assign PDB residues to UniProt/Pfam
positions.  
Given the fact that pairwise alignments can lead to misleading association between
residues in both sequences, SIFTS offers  more reliable association between sequence and
structure residue numbers.  

```julia
using MIToS.SIFTS # to load the SIFTS module
```  

## Features

- Download and parse SIFTS XML files
- Store residue-level mapping in Julia
- Easy generation of `Dict`s between residues numbers

## Contents

```@contents
Pages = ["SIFTS.md"]
Depth = 4
```  

## [Simplest residue-level mapping](@id Simplest-residue-level-mapping)

This module export the function `siftsmapping` to generate a `Dict` between residue
numbers. This function takes 5 positional arguments.
    1) The name of the SIFTS XML file to parse,  
    2) the source database,  
    3) the source protein/structure identifier,  
    4) the destiny database and,  
    5) the destiny protein/structure identifier.  
Optionally it’s possible to indicate a particular PDB `chain` and if `missings` will be used.  

Databases should be indicated using an available sub-type of `DataBase`. Keys and values
types will be depend on the residue number type in that database.

| Type `db...`  | Database | Residue number type |
|---------------|----------|---------------------|
| `dbPDBe`		| **PDBe** (Protein Data Bank in Europe) | `Int` |
| `dbInterPro`	| **InterPro** | `ASCIIString` |
| `dbUniProt`	| **UniProt** | `Int` |
| `dbPfam`		| **Pfam** (Protein families database) | `Int` |
| `dbNCBI`		| **NCBI** (National Center for Biotechnology Information) | `Int` |
| `dbPDB`		| **PDB** (Protein Data Bank) | `ASCIIString` |
| `dbCATH`		| **CATH** | `ASCIIString` |
| `dbSCOP` 		| **SCOP** (Structural Classification of Proteins) | `ASCIIString` |

To download the XML SIFTS file of a determined PDB use the `downloadsifts` function.  

```@example sifts_simple
using MIToS.SIFTS

siftsfile = downloadsifts("1IVO")
```  

The following example, shows the residue number mapping between *Pfam* and *PDB*.
*Pfam* uses *UniProt* coordinates and *PDB* uses their own residue numbers with insertion
codes. Note that **the `siftsmapping` function is case sensitive**, and that
**SIFTS stores PDB identifiers using lowercase characters**.  


```@example sifts_simple
siftsmap = siftsmapping(siftsfile,
                        dbPfam,
                        "PF00757",
                        dbPDB,
                        "1ivo", # SIFTS stores PDB identifiers in lowercase
                        chain="A", # In this example we are only using the chain A of the PDB
                        missings=false) # Residues without coordinates aren't used in the mapping
```  

## [Storing residue-level mapping](@id Storing-residue-level-mapping)

If you need more than the residue number mapping between two databases, you could access
all the residue-level cross references using the function `read` in the `SIFTSXML``Format`
file. The `parse` function (and therefore the `read` function) for the `SIFTSXML` format,
also takes the keyword arguments `chain` and `missings`. The `read`/`parse` function
returns a `Vector` of `SIFTSResidue`s objects that stores the cross references between
residues in each database.  

```@example sifts_simple
siftsresidues = read(siftsfile, SIFTSXML, chain="A", missings=false) # Array{SIFTSResidue,1}

residue_data = siftsresidues[300]
```

You are free to access the `SIFTSResidue` fields in order to get the desired information.
`SIFTSResidue` objects contain `db...` objects (sub-types of `DataBase`), with the cross
referenced information. You should note that, except for the `PDBe` and `InterPro` fields,
the fields are `Nullable`s objects so, you need to use the `get` function to access the
`db...` object. For example, getting the UniProt residue name
(one letter code of the amino acid) would be:  

```@example sifts_simple
isnull(residue_data.UniProt) ? "" : get(residue_data.UniProt).name
```  

That line of code returns an empty string if the UniProt field is null. Otherwise, it
returns a string with the name of the residue in UniProt. Because that way of access
values in a Residue is too verbose, MIToS defines a more complex signature for `get`.
Using MIToS `get` the previous line of code will be:  

```@example sifts_simple
#   SIFTSResidue  database   field  default
get(residue_data, dbUniProt, :name, "")
```  

The is not need to use the full signature, but the returned value will change. In
particular, a `Nullable` object is returned if a default value is not given at the end of
the signature:  

```@repl
using MIToS.SIFTS; residue_data = read("1ivo.xml.gz", SIFTSXML)[301] # hide
get(residue_data, dbUniProt) # Takes the database type and returns a nullable with the field content
get(residue_data, dbUniProt, :name) # Takes also a Symbol with a field name and returns a nullable with the field content inside the database type
```  

But you don't need the `get`function to access the three letter code of the residue in
`PDBe` because the `PDBe` field is not `Nullable`.

```@example sifts_simple
residue_data.PDBe.name
```  

`SIFTSResidue` also store information about if that residue is `missing` in the
PDB structure and the information about the secondary structure (`sscode` and `ssname`):  

```@repl
using MIToS.SIFTS; residue_data = read("1ivo.xml.gz", SIFTSXML)[301] # hide
residue_data.missing
residue_data.sscode
residue_data.ssname
```  

### [Accessing residue-level cross references](@id Accessing-residue-level-cross-references)

You can ask for particular values in a single `SIFTSResidue` using the `get` function.  

```@repl
using MIToS.SIFTS
residue_data = read("1ivo.xml.gz", SIFTSXML)[301]
# Is the UniProt residue name in the list of basic amino acids ["H", "K", "R"]?
get(residue_data, dbUniProt, :name, "") in ["H", "K", "R"]
```  

Use higher order functions and lambda expressions (anonymous functions) or
list comprehension to easily ask for information on the `Vector{SIFTSResidue}`. You can
use `get` with the previous signature or simple `get`, direct field access and `isnull`.

```@example sifts_simple
# Captures PDB residue numbers if the Pfam id is "PF00757"
resnums = [ get(res.PDB).number for res in siftsresidues if !isnull(res.PDB) && get(res, dbPfam, :id, "") == "PF00757" ]
```  

**Useful higher order functions are:**

**`find`**  

```@example sifts_simple
# Which of the residues have UniProt residue names in the list ["H", "K", "R"]? (basic residues)
indexes = find(res -> get(res, dbUniProt, :name, "") in ["H", "K", "R"], siftsresidues)
```  

**`map`**  

```@example sifts_simple
map(i -> get(siftsresidues[i].UniProt), indexes) # UniProt data of the basic residues
```

**`filter`**  

```@example sifts_simple
# SIFTSResidues with UniProt names in ["H", "K", "R"]
basicresidues = filter(res -> get(res, dbUniProt, :name, "") in ["H", "K", "R"], siftsresidues)

get(basicresidues[1].UniProt) # UniProt data of the first basic residue
```  

#### [Example: Which residues are missing in the PDB structure](@id Example:-Which-residues-are-missing-in-the-PDB-structure)

Given that `SIFTSResidue` objects store a `missing` residue flag, it’s easy to get a
vector where there is a `true` value if the residue is missing in the structure.  

```@repl
using MIToS.SIFTS
sifts_1ivo = read("1ivo.xml.gz", SIFTSXML, chain="A"); # SIFTSResidues of the 1IVO chain A
[res.missing for res in sifts_1ivo]
```  

However, if you need to filter using other conditions, you’ll find useful the `get`
function. In this example, we are going to ask for the *UniProt id*
(to avoid problems with fragments, tags or chimeric/fusion proteins). We are also using
`get` to select an specific PDB chain.  

```@setup sifts_1jqz
using MIToS.SIFTS
downloadsifts("1JQZ")
```

```@repl
using MIToS.SIFTS
sifts_1jqz = read("1jqz.xml.gz", SIFTSXML); # It has an amino terminal his tag
missings = [ (  ( get(res, dbUniProt, :id, "") == "P05230" ) &
                ( get(res, dbPDB, :chain, "") ==  "A" ) &
                res.missing ) for res in sifts_1jqz             ];
println("There are only ", sum(missings), " missing residues in the chain A, associated to UniProt P05230")
println("But there are ", sum([ res.missing for res in sifts_1jqz ]), " missing residues in the PDB file.")
```  
