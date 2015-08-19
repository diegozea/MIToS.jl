module SIFTS

  using LightXML
  using MIToS.Utils

  export downloadsifts, siftsmapping, siftsresidues

  include("XMLParser.jl")
  include("ResidueMapping.jl")

end
