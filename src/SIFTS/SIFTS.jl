"""
The `SIFTS` module of MIToS allows to obtain the
residue-level mapping between databases stored in the SIFTS XML files.
It makes easy to assign PDB residues to UniProt/Pfam positions.
Given the fact that pairwise alignments can lead to misleading association between residues in both sequences,
SIFTS offers  more reliable association between sequence and structure residue numbers.

**Features**

  - Download and parse SIFTS XML files
  - Store residue-level mapping in Julia
  - Easy generation of `OrderedDict`s between residues numbers

```julia
using MIToS.SIFTS
```
"""
module SIFTS

import LightXML

using AutoHashEquals
using OrderedCollections
using MIToS.Utils

export DataBase,
    dbPDBe,
    dbInterPro,
    dbUniProt,
    dbPfam,
    dbNCBI,
    dbPDB,
    dbCATH,
    dbSCOP,
    dbSCOP2,
    dbSCOP2B,
    dbEnsembl,
    SIFTSResidue,
    downloadsifts,
    siftsmapping,
    SIFTSXML,
    # Mitos.Utils
    All,
    read_file,
    parse_file,
    # Imported from Base (and exported for docs)
    parse

include("XMLParser.jl")
include("ResidueMapping.jl")

end
