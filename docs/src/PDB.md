```@setup log
@info "PDB docs"
```

# [PDB](@id Module-PDB)

The module `PDB` defines types and methods to work with protein structures inside Julia. It
is useful to link structural and sequential information, and needed for measure the
predictive performance at protein contact prediction of mutual information scores.  

```julia
using MIToS.PDB # to load the PDB module
```  

## Features  

- [**Read and parse**](@ref Read-and-parse-PDB-files) PDB and PDBML files.
- Calculate distance and contacts between atoms or residues.
- Determine interaction between residues.

## Contents  

```@contents
Pages = ["PDB.md"]
Depth = 4
```  

## Retrieve information from PDB database  

This module exports the `downloadpdb` function, to retrieve a PDB file from  
[PDB database![](./assets/external-link.png)](http://www.rcsb.org/pdb/home/home.do). This
function downloads a gzipped `PDBML` file, which could be easily read it with MIToS
by default, but you are able to determine the `format` as `PDBFile` if you want it.  

```@example pdb_io
using MIToS.PDB

pdbfile = downloadpdb("1IVO", format=PDBFile)
```  

`PDB` module also exports a `getpdbdescription` to access the header information of a
PDB entry.  

```@example pdb_io
getpdbdescription("1IVO")
```  

## [Read and parse PDB files](@id Read-and-parse-PDB-files)  

This is easy using the `read` and `parse` functions, indicating the filename and the
`FileFormat`: `PDBML` for PDB XML files or `PDBFile` for usual PDB files. These functions
returns a `Vector` of `PDBResidue` objects with all the residues in the PDB.  
To return only a specific subset of residues/atoms you can use any of the following
keyword arguments:  

|keyword arguments | default | returns only ... |
|------------------|---------|-----------|
|`chain` | `All` | residues from a PDB chain, i.e. `"A"` |
|`model` | `All` | residues from a determined model, i.e. `"1"` |
|`group` | `All` | residues from a group: `"ATOM"`, `"HETATM"` or `All` for both |
|`atomname` | `All` | atoms with a specific name, i.e. `"CA"` |
|`onlyheavy` | `false` | heavy atoms (not hydrogens) if it's `true` |
|`occupancyfilter` | `false` | only the atoms with the best occupancy are returned if it's `true` |

!!! note
    **For PDBML files** it is possible to use the keyword argument `label` to `false`
    (default to `true`) to get the **auth_** attributes instead of the **label_**
    attributes for `chain`, `atom` and residue `name` fields. The **auth_** attributes are
    alternatives provided by an author in order to match the identification/values used
    in the publication that describes the structure.  

```@example pdb_io
# Read α carbon of each residue from the 1ivo pdb file, in the model 1, chain A and in the ATOM group.
CA_1ivo = read(pdbfile, PDBFile, model="1", chain="A", group="ATOM", atomname="CA")

CA_1ivo[1] # First residue. It has only the α carbon.
```

## Looking for particular residues  

MIToS parse PDB files to vector of residues, instead of using a hierarchical structure
like other packages. This approach makes the search and selection of residues or atoms a
little different.
To make it easy, this module exports a number of functions and macros to select particular
residues or atoms. Given the fact that residue numbers from different chains, models, etc.
can collide, **it's mandatory to indicate the `model`, `chain`, `group`, `residue` number
and `atom` name in a explicit way** to these functions or macros. If you want to select all
the residues in one of the categories, you are able to use the type `All`. You can also use
regular expressions or functions to make the selections.

```@example pdb_select
using MIToS.PDB
pdbfile = downloadpdb("1IVO", format=PDBFile)
residues_1ivo = read(pdbfile, PDBFile)
# Select residue number 9 from model 1 and chain B
residues(residues_1ivo, "1", "B", All, "9")
```

### Getting a `Dict` of `PDBResidue`s

If you prefer a `Dict` of `PDBResidue`, indexed by their residue numbers, you can use the
`residuedict` function or the `@residuedict` macro.  

```@example pdb_select
# Dict of residues from the model 1, chain A and from the ATOM group
chain_a = residuesdict(residues_1ivo, "1", "A", "ATOM", All)
chain_a["9"]
```  

You can do the same with the macro `@residuesdict` to get a more readable code  

```@example pdb_select
chain_a = @residuesdict residues_1ivo model "1" chain "A" group "ATOM" residue All
chain_a["9"]
```  

### Select particular residues  

Use the `residues` function to collect specific residues. It's possible to use a single
**residue number** (i.e. `"2"`) or even a **function** which should return true for the
selected residue numbers. Also **regular expressions** can be used to select residues.
Use `All` to select all the residues.  


```@example pdb_select
residue_list = map(string, 2:5)

# If the list is large, you can use a `Set` to gain performance
# residue_set = Set(map(string, 2:5))
```

```@example pdb_select
first_res = residues(residues_1ivo, "1", "A", "ATOM", resnum -> resnum in residue_list)

for res in first_res
    println(res.id.name, " ", res.id.number)
end
```

A more complex example using an anonymous function:  

```@example pdb_select
# Select all the residues of the model 1, chain A of the ATOM group with residue number less than 5

first_res = residues(residues_1ivo, "1", "A", "ATOM", x -> parse(Int, match(r"^(\d+)", x)[1]) <= 5 )
# The anonymous function takes the residue number (string) and use a regular expression
# to extract the number (without insertion code).
# It converts the number to `Int` to test if the it is `<= 5`.

for res in first_res
    println(res.id.name, " ", res.id.number)
end
```

Use the `@residues` macro for a cleaner syntax.  

```@example pdb_select
# You can use All, regular expressions or functions also for model, chain and group:

# i.e. Takes the residue 10 from chains A and B

for res in @residues residues_1ivo model "1" chain ch -> ch in ["A","B"] group "ATOM" residue "10"
    println(res.id.chain, " ", res.id.name, " ", res.id.number)
end
```

### Select particular atoms

The `atoms` function or macro allow to select a particular set of atoms.

```@example pdb_select
# Select all the atoms with name starting with "C" using a regular expression
# from all the residues of the model 1, chain A of the ATOM group

carbons = @atoms residues_1ivo model "1" chain "A" group "ATOM" residue All atom r"C.+"

carbons[1]
```  

You can also use the `atoms` function instead of the `@atoms` macro:  

```@example pdb_select
atoms(residues_1ivo, "1", "A", "ATOM", All, r"C.+")[1]
```

## Protein contact map

The PDB module offers a number of functions to measure `distance`s between atoms or
residues, to detect possible interactions or `contact`s. In particular the `contact`
function calls the `distance` function using a threshold or limit in an optimized way.
The measure can be done between alpha carbons (`"CA"`), beta carbons (`"CB"`) (alpha carbon
for glycine), any heavy atom (`"Heavy"`) or any (`"All"`) atom of the residues.

In the following **example**, whe are going to plot a contact map for the *1ivo* chain A.
Two residues will be considered in contact if their β carbons (α carbon for glycine) have a
distance of 8Å or less.  

```@example pdb_cmap
using MIToS.PDB

pdbfile = downloadpdb("1IVO", format=PDBFile)

residues_1ivo = read(pdbfile, PDBFile)

pdb = @residues residues_1ivo model "1" chain "A" group "ATOM" residue All

dmap = distance(pdb, criteria="All") # Minimum distance between residues using all their atoms
```

Use the `contact` function to get a contact map:  

```@example pdb_cmap
cmap = contact(pdb, 8.0, criteria="CB") # Contact map
```

```@setup pdb_cmap
@info "PDB: Cmap"
using Plots
gr() # Hide possible warnings
```

```@example pdb_cmap
using Plots
gr()

heatmap(dmap, grid=false, yflip=true, ratio=:equal)

png("pdb_dmap.png") # hide
nothing # hide
```  

![](pdb_dmap.png)  


```@example pdb_cmap
heatmap(cmap, grid=false, yflip=true, ratio=:equal)

png("pdb_cmap.png") # hide
nothing # hide
```  

![](pdb_cmap.png)  

## Structural superposition  

```@setup pdb_rmsd
@info "PDB: RMSD"
using Plots
gr() # Hide possible warnings
```

```@example pdb_rmsd
using MIToS.PDB

pdbfile = downloadpdb("2HHB")

res_2hhb = read(pdbfile, PDBML)

chain_A = pdb = @residues res_2hhb model "1" chain "A" group "ATOM" residue All
chain_C = pdb = @residues res_2hhb model "1" chain "C" group "ATOM" residue All

using Plots
gr()

scatter3d(chain_A, label="A", alpha=0.5)
scatter3d!(chain_C, label="C", alpha=0.5)

png("pdb_unaligned.png") # hide
nothing # hide
```  

![](pdb_unaligned.png)  

```@example pdb_rmsd
superimposed_A, superimposed_C, RMSD = superimpose(chain_A, chain_C)

RMSD
```

```@example pdb_rmsd
scatter3d(superimposed_A, label="A", alpha=0.5)
scatter3d!(superimposed_C, label="C", alpha=0.5)
png("pdb_aligned.png") # hide
nothing # hide
```  

![](pdb_aligned.png)  
