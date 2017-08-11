
# SIFTS 

The `SIFTS` module of MIToS allows to obtain the residue-level mapping between databases stored in the SIFTS XML files. It makes easy to assign PDB residues to UniProt/Pfam positions.  
Given the fact that pairwise alignments can lead to misleading association between residues in both sequences, SIFTS offers  more reliable association between sequence and structure residue numbers.


## Features

- Download and parse SIFTS XML files
- Store residue-level mapping in Julia
- Easy generation of `Dict`s between residues numbers



```julia
using MIToS.SIFTS
```

## Contents

- [Simplest residue-level mapping](#Simplest-residue-level-mapping)
- [Storing residue-level mapping](#Storing-residue-level-mapping)
- [Accessing residue-level cross references](#Accessing-residue-level-cross-references)
    - [capture](#capture) 
    - [collectcaptures](#collectcaptures) 
    - [isobject](#isobject)
    - [findobjects](#findobjects)
    - [collectobjects](#collectobjects)
- [Example: Which residues are missing in the PDB structure](#Example:-Which-residues-are-missing-in-the-PDB-structure)


```julia
# Truncate IJulia outputs at:
ENV["LINES"]   = 15 
ENV["COLUMNS"] = 60;
```

## Simplest residue-level mapping  

This module export the function `siftsmapping` to generate a `Dict` between residue numbers. This function takes 5 positional arguments. 1) The name of the SIFTS XML file to parse, 2) the source database 3) the source protein/structure identifier, 4) the destiny database and 5) the destiny protein/structure identifier. Optionally it’s possible to indicate a particular PDB `chain` and if `missings` will be used.  

Databases should be indicated using an available sub-type of `DataBase`. Keys and values types will be depend on the residue number type in that database.

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


```julia
siftsfile = downloadsifts("1IVO")
```

      % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                     Dload  Upload   Total   Spent    Left  Speed
    100 53962  100 53962    0     0   7749      0  0:00:06  0:00:06 --:--:-- 11338





    "1ivo.xml.gz"



The following example, shows the residue number mapping between *Pfam* and *PDB*. *Pfam* uses *UniProt* coordinates and *PDB* uses their own residue numbers with insertion codes. Note that <span class="text-warning">the `siftsmapping` function is case sensitive</span>, and that <span class="text-warning">SIFTS stores PDB identifiers using lowercase characters</span>.  


```julia
siftsmap = siftsmapping( siftsfile, 
                    dbPfam, "PF00757", 
                    dbPDB, "1ivo", # SIFTS stores PDB identifiers in lowercase
                    chain="A", # In this example we are only using the chain A of the PDB
                    missings=false) # Residues without coordinates aren't used in the mapping
```




    Dict{Int64,ASCIIString} with 162 entries:
      329 => "305"
      210 => "186"
      288 => "264"
      241 => "217"
      267 => "243"
      306 => "282"
      275 => "251"
      197 => "173"
      215 => "191"
      181 => "157"
      ⋮   => ⋮



<a href="#"><i class="fa fa-arrow-up"></i></a>

## Storing residue-level mapping  

If you need more than the residue number mapping between two databases, you could access all the residue-level cross references using the function `read` in the `SIFTSXML``Format` file. The `parse` function (and therefore the `read` function) for the `SIFTSXML` format, also takes the keyword arguments `chain` and `missings`. The `read`/`parse` function returns a `Vector` of `SIFTSResidue`s objects that stores the cross references between residues in each database.  


```julia
siftsresidues = read(siftsfile, SIFTSXML, chain="A", missings=false) # Array{SIFTSResidue,1}

residue_data = siftsresidues[300]
```




    SIFTSResidue
      PDBe:
        number: 301
        name: LYS
      UniProt (Nullable) :
        id: P00533
        number: 325
        name: K
      Pfam (Nullable) :
        id: PF00757
        number: 325
        name: K
      NCBI (Nullable) :
        id: 9606
        number: 325
        name: K
      PDB (Nullable) :
        id: 1ivo
        number: 301
        name: LYS
        chain: A
      SCOP (Nullable) :
        id: 76847
        number: 301
        name: LYS
        chain: A
      CATH (Nullable) :
        id: 2.10.220.10
        number: 301
        name: LYS
        chain: A
        InterPro: [MIToS.SIFTS.dbInterPro("IPR009030","325","K","SSF57184")]




You are free to access the `SIFTSResidue` fields in order to get the desired information. `SIFTSResidue` objects contain `db...` objects (sub-types of `DataBase`), with the cross referenced information. You should note that, except for the `PDBe` and `InterPro` fields, the fields are `Nullable`s objects so, you need to use the `get` function to access the `db...` object. For example, getting the UniProt residue name (one letter code of the amino acid) would be:


```julia
get(residue_data.UniProt).name
```




    "K"



But you don't need the `get`function to access the three letter code of the residue in `PDBe`


```julia
residue_data.PDBe.name
```




    "LYS"



`SIFTSResidue` also store information about if that residue is `missing` in the PDB structure:


```julia
residue_data.missing
```




    false



<a href="#"><i class="fa fa-arrow-up"></i></a>

### Accessing residue-level cross references

To easily ask for information to the `Vector{SIFTSResidue}` use the following functions and the tests: `Is`, `In` and/or `Not`.

#### `capture`

Takes a `SIFTSResidue`, a `db...` type and a `Symbol` with the name of the field to capture from that database. Returns a `Nullable` with the field content if the test are passed over a determined database.


```julia
# Captures the residue name in UniProt if the residue number in PDB is "301"
capture(residue_data, dbUniProt, :name, dbPDB, Is(:number, "301"))
```




    Nullable("K")



#### `collectcaptures`

Returns a vector of `Nullable`s with the `capture`s of the `field`s. The element is null if any test fails or the object hasn't the `field`.


```julia
# Captures PDB residue numbers if the Pfam id is "PF00757"
captures = collectcaptures(siftsresidues, dbPDB, :number, dbPfam, Is(:id, "PF00757"))
```




    511-element Array{Nullable{ASCIIString},1}:
     Nullable{ASCIIString}()
     Nullable{ASCIIString}()
     Nullable{ASCIIString}()
     Nullable{ASCIIString}()
     Nullable{ASCIIString}()
     ⋮                      
     Nullable{ASCIIString}()
     Nullable{ASCIIString}()
     Nullable{ASCIIString}()
     Nullable{ASCIIString}()
     Nullable{ASCIIString}()




```julia
captures[ [!isnull(res) for res in captures] ] # Selects not null elements
```




    162-element Array{Nullable{ASCIIString},1}:
     Nullable("153")
     Nullable("154")
     Nullable("155")
     Nullable("156")
     Nullable("157")
     ⋮              
     Nullable("310")
     Nullable("311")
     Nullable("312")
     Nullable("313")
     Nullable("314")



#### `isobject` 

Returns `true` if the tests are successfully passed for that `DataBase` sub-type on that `SIFTSResidue`.


```julia
# Is it a basic amino acid?
# Is the UniProt residue name in the list ["H", "K", "R"]
isobject(residue_data, dbUniProt, SIFTS.In(:name, ["H", "K", "R"]))
```




    true



#### `findobjects`

Returns a vector of the indexes for which `isobject` is `true` in the input vector.


```julia
# Which of the residues are basic?
# Which of the residues have UniProt residue names in the list ["H", "K", "R"]?
indexes = findobjects(siftsresidues, dbUniProt, SIFTS.In(:name, ["H", "K", "R"]))
```




    69-element Array{Int64,1}:
       3
       4
      12
      22
      28
       ⋮
     482
     496
     502
     506
     508




```julia
[ get(siftsresidues[ idx ].UniProt) for idx in indexes ] # UniProt data of the basic residues
```




    69-element Array{Any,1}:
     MIToS.SIFTS.dbUniProt("P00533",28,"K") 
     MIToS.SIFTS.dbUniProt("P00533",29,"K") 
     MIToS.SIFTS.dbUniProt("P00533",37,"K") 
     MIToS.SIFTS.dbUniProt("P00533",47,"H") 
     MIToS.SIFTS.dbUniProt("P00533",53,"R") 
     ⋮                                      
     MIToS.SIFTS.dbUniProt("P00533",507,"H")
     MIToS.SIFTS.dbUniProt("P00533",521,"R")
     MIToS.SIFTS.dbUniProt("P00533",527,"R")
     MIToS.SIFTS.dbUniProt("P00533",531,"R")
     MIToS.SIFTS.dbUniProt("P00533",533,"R")



#### `collectobjects`

Returns a vector with the objects for which `isobject` is `true`.


```julia
# Collect SIFTSResidues with UniProt names in ["H", "K", "R"]
basicresidues = collectobjects(siftsresidues, dbUniProt, SIFTS.In(:name, ["H", "K", "R"]))

get(basicresidues[1].UniProt) # UniProt data of the first basic residue
```




    MIToS.SIFTS.dbUniProt("P00533",28,"K")



<a href="#"><i class="fa fa-arrow-up"></i></a>

#### Example: Which residues are missing in the PDB structure

Given that `SIFTSResidue` objects store a `missing` residue flag, it’s easy to get a vector where there is a `true` value if the residue is missing in the structure.


```julia
sifts_1ivo = read(siftsfile, SIFTSXML, chain="A") # Array{SIFTSResidue,1} of the 1IVO chain A

[res.missing for res in sifts_1ivo]
```




    622-element Array{Any,1}:
      true
     false
     false
     false
     false
         ⋮
      true
      true
      true
      true
      true



However, if you need to filter using other conditions, you’ll find useful the `isobject` function. In this example, we are going to ask for the *UniProt id* (to avoid problems with fragments, tags or chimeric/fusion proteins). We are also using `isobject` to select an specific PDB chain.


```julia
file_1jqz = downloadsifts("1JQZ")
```

      % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                     Dload  Upload   Total   Spent    Left  Speed
    100 21058  100 21058    0     0   2611      0  0:00:08  0:00:08 --:--:--  5709





    "1jqz.xml.gz"




```julia
sifts_1jqz = read(file_1jqz, SIFTSXML) # It has an amino terminal his tag

missings = [ isobject(res, dbUniProt, Is(:id, "P05230")) & 
             isobject(res, dbPDB, Is(:chain, "A")) & 
             res.missing for res in sifts_1jqz             ]

println("There are only ", sum(missings), " missing residues in the chain A, associated to UniProt P05230")
```

    There are only 3 missing residues in the chain A, associated to UniProt P05230



```julia
println("But there are ", sum([ res.missing for res in sifts_1jqz ]), " missing residues in the PDB file.")
```

    But there are 10 missing residues in the PDB file.

