"""
The `Utils` has common utils functions and types used in other modules.

```julia
using MIToS.Utils
```
"""
module Utils

using Downloads
using CodecZlib
using NamedArrays
using Logging

import LightXML

export  # GeneralUtils.jl
    All,
    get_n_words,
    hascoordinates,
    select_element,
    matrix2list,
    list2matrix,
    check_pdbcode,
    getarray,
    # Read.jl
    FileFormat,
    lineiterator,
    check_file,
    isnotemptyfile,
    download_file,
    read_file,
    parse_file,
    # Write.jl
    Commandline,
    write_file,
    print_file,
    # ThreeLetterResidues.jl
    THREE2ONE,
    # Imported from Base (and exported for docs)
    read,
    write


include("GeneralUtils.jl")
include("Read.jl")
include("Write.jl")
include("ThreeLetterResidues.jl")

@deprecate deleteitems!(vector::Vector, items) filter!(x -> x ∉ items, vector)

# COMMAND LINE: Scripts Module
# ============================

module Scripts

    export run_single_script,
        script,
        open_output,
        close_output,
        readorparse,
        loadedversion

    using Pkg
    using CodecZlib
    using MIToS.Utils # to use read_file

    include("Scripts.jl")

end # Scripts

end # Utils
