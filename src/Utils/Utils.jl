module Utils

  export deleteitems!, get_n_words, hascoordinates, select_element, matrix2list,

  # eachline,

  Format,

  AbstractTest, Is, In, Not, capture, isobject, findobjects, collectobjects, collectcaptures, guess_type

  include("generalutils.jl")
  include("EachLineString.jl")
  include("Read.jl")
  include("FindObjects.jl")

end
