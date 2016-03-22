module Utils

export deleteitems!, get_n_words, hascoordinates, select_element, matrix2list, list2matrix, check_file, isnotemptyfile,

# eachline,

Format,

AbstractTest, TestType, TestOperation, Is, In, Not, capture, isobject, findobjects, collectobjects, collectcaptures, guess_type,

Commandline

include("generalutils.jl")
include("EachLineString.jl")
include("Read.jl")
include("Write.jl")
include("FindObjects.jl")


# COMMAND LINE
# ============

module Commandline

export parse_commandline, runfun, set_parallel, open_output, close_output

using ArgParse

include("Commandline.jl")

end # Commandline

end # Utils
