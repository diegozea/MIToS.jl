module SIFTS

  using LightXML
  using AutoHashEquals
  using MIToS.Utils

  import Base: hash, ==, call, string, print, write, show, convert, isnull
  import MIToS.Utils: isobject, findobjects, collectobjects, capture, collectcaptures, guess_type

  export DataBase, dbPDBe, dbInterPro, dbUniProt, dbPfam, dbNCBI, dbPDB, dbCATH, dbSCOP,
  SIFTSResidue, has, getdatabase, getcoordinate, ischain,
  downloadsifts, siftsmapping, siftsresidues,

  # Mitos.Utils
  capture, collectcaptures, isobject, findobjects, collectobjects, Is, Not, In

  include("XMLParser.jl")
  include("ResidueMapping.jl")

end
