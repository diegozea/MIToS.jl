#!/usr/bin/env julia

using MIToS.Utils.Commandline

Args = parse_commandline(
# TO DO -----------------------------------------------------------------------
    ["--arg", "-a"],
    Dict(
        :help => "Argument",
        :arg_type => Int,
        :default => 0
    ),
    # Keywords...
    description="Made with MIToS",
    output=".mitos."
# -----------------------------------------------------------------------------
    )

set_parallel(Args["parallel"])

@everywhere begin

    const args = remotecall_fetch(1,()->Args)
    const file = args["FILE"]
    const list = args["list"]
    const output = args["output"]

    # TO DO -----------------------------------------------------------------------
    using MIToS
    # -----------------------------------------------------------------------------

    function main(input) # input must be a file or STDIN from FILE (call_main)
        fh_out = open_output(input, Main.output)
        try
            # TO DO -----------------------------------------------------------------------
            arg_one = args["arg"]
            println(fh_out, "RUN : $arg_one : $input")
            # -----------------------------------------------------------------------------
        catch err
            println("ERROR: ", input)
            println(err)
        finally
            close(fh_out)
        end
    end

end

runfun(main, file, list)
