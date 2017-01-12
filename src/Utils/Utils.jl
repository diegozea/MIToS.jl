module Utils

using Requests
using GZip
using LightXML

export  # GeneralUtils.jl
        All,
        get_n_words,
        hascoordinates,
        select_element,
        matrix2list, list2matrix,
        check_pdbcode,
        # Read.jl
        Format,
        lineiterator,
        check_file, isnotemptyfile,
        download_file,
        # Write.jl
        Commandline

include("GeneralUtils.jl")
include("Read.jl")
include("Write.jl")

@deprecate deleteitems!(vector::Vector, items) filter!(x -> x ∉ items, vector)

# COMMAND LINE: Scripts Module
# ============================

module Scripts

export  parse_commandline,
        runscript, run_single_script,
        script,
        set_parallel,
        open_output, close_output,
        readorparse

using ArgParse, GZip

include("Scripts.jl")

end # Scripts

end # Utils
