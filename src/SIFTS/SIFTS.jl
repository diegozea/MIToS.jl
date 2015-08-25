module SIFTS

  using LightXML
  using AutoHashEquals
  using MIToS.Utils

  import Base: hash, ==, call, string, print, write, show, convert, isnull

  export DataBase, dbPDBe, dbInterPro, dbUniProt, dbPfam, dbNCBI, dbPDB, dbCATH, dbSCOP,
  SIFTSResidue, has, getdatabase, getcoordinate, ischain,
  downloadsifts, siftsmapping, siftsresidues

  include("XMLParser.jl")
  include("ResidueMapping.jl")

end
