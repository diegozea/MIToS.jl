module Utils

  export deleteitems!, get_n_words, hascoordinates, select_element, matrix2list, list2matrix,

  # eachline,

  Format,

  AbstractTest, TestType, TestOperation, Is, In, Not, capture, isobject, findobjects, collectobjects, collectcaptures, guess_type

  include("generalutils.jl")
  include("EachLineString.jl")
  include("Read.jl")
  include("FindObjects.jl")

end
