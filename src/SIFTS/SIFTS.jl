module SIFTS

  using LightXML
  using AutoHashEquals
  using MIToS.Utils

  import Base: hash, ==

  export CoordinateSystem,
  PDBeCoordinate, UniProtCoordinate, PDBresnumCoordinate,
  DataBase,
  RefPDBe, RefUniProt, RefPfam, RefNCBI, RefPDB, RefCATH, RefSCOP,
  dbPDBe, dbUniProt, dbPfam, dbNCBI, dbPDB, dbCATH, dbSCOP,

  SIFTSResidue, has, getdatabase, getcoordinate, ischain,

  downloadsifts, siftsmapping, siftsresidues

  include("XMLParser.jl")
  include("ResidueMapping.jl")

end
