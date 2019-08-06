#!/usr/bin/env julia

using Distributed
using MIToS.Utils.Scripts

Args = parse_commandline(
# TO DO -----------------------------------------------------------------------
    description="""
    Creates a file in Stockholm format with the aligned columns from a Pfam Stockholm file.
    Insertions are deleted, as they are unaligned in a proﬁle HMM.
    The output file *.aligned.* contains UniProt residue numbers and original column numbers in its annotations.
    """,
    output=".aligned."
# -----------------------------------------------------------------------------
    )

set_parallel(Args["parallel"])

@everywhere begin

    const args = remotecall_fetch(()->Args,1)

    import MIToS.Utils.Scripts: script

    # TO DO -----------------------------------------------------------------------
    using MIToS.MSA

    function script(input::Union{Base.LibuvStream,  AbstractString},
                    args,
                    fh_out::Union{Base.LibuvStream, IO})
        try
            aln = readorparse(input, Stockholm, generatemapping=true, useidcoordinates=true, deletefullgaps=true)
            print(fh_out, aln, Stockholm)
        catch err
            @warn(string("ERROR for ", input, ": ", err))
        end

    end
    # -----------------------------------------------------------------------------

end

runscript(args)
