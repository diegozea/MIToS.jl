#!/usr/bin/env julia

using Distributed
using MIToS.Utils.Scripts

Args = parse_commandline(
    # TO DO ----------------------------------------------------------------------
    ["--arg", "-a"],
    Dict(
        :help => "Argument",
        :arg_type => Int,
        :default => 0
    ),
    # Keywords...
    description="Made with MIToS",
    output=".mitos."
    # ----------------------------------------------------------------------------
    )

set_parallel(Args["parallel"])

@everywhere begin

    const args = remotecall_fetch(()->Args,1)

    import MIToS.Utils.Scripts: script

    # TO DO ----------------------------------------------------------------------
    using MIToS
    # ----------------------------------------------------------------------------

    function script(input::Union{Base.LibuvStream,  AbstractString},
                    args,
                    fh_out::Union{Base.LibuvStream, IO})
        # TO DO ------------------------------------------------------------------
        arg_one = args["arg"]
        println(fh_out, "RUN : $arg_one")
        println("$input : ")
        println(read(input, String))
        # ------------------------------------------------------------------------
    end

end

runscript(args)
