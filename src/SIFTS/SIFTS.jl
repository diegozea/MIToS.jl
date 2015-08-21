module SIFTS

  using LightXML
  using MIToS.Utils

  export CoordinateSystem,
  PDBeCoordinate, UniProtCoordinate, PDBresnumCoordinate,
  DataBase,
  dbPDBe, dbUniProt, dbPfam, dbInterPro, dbNCBI, dbPDB, dbCATH, dbSCOP,

  downloadsifts, siftsmapping, siftsresidues

  include("XMLParser.jl")
  include("ResidueMapping.jl")

end
