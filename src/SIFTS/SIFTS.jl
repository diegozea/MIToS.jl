"""
The `SIFTS` module of MIToS allows to obtain the
residue-level mapping between databases stored in the SIFTS XML files.
It makes easy to assign PDB residues to UniProt/Pfam positions.
Given the fact that pairwise alignments can lead to misleading association between residues in both sequences,
SIFTS offers  more reliable association between sequence and structure residue numbers.


**Features**

- Download and parse SIFTS XML files
- Store residue-level mapping in Julia
- Easy generation of `Dict`s between residues numbers

```julia

using MIToS.SIFTS
```
"""
module SIFTS

  using LightXML
  using AutoHashEquals
  using MIToS.Utils

  import Base: hash, ==, call, string, print, write, show, convert, isnull, parse
  import MIToS.Utils: isobject, findobjects, collectobjects, capture, collectcaptures, guess_type

  export DataBase, dbPDBe, dbInterPro, dbUniProt, dbPfam, dbNCBI, dbPDB, dbCATH, dbSCOP,
  SIFTSResidue, downloadsifts, siftsmapping, SIFTSXML,

  # Mitos.Utils
  capture, collectcaptures, isobject, findobjects, collectobjects, Is, Not, In

  include("XMLParser.jl")
  include("ResidueMapping.jl")

end
