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


# COMMAND LINE: Scripts Module
# ============================

module Scripts

export parse_commandline,
runscript, run_single_script, script,
set_parallel, open_output, close_output,
readorparse

using ArgParse, GZip

include("scripts.jl")

end # Scripts

end # Utils
