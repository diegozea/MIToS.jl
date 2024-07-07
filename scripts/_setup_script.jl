# Environment setup
# -----------------
using Pkg

# MIToS is imported before including this file
const SCRIPTS_PATH = joinpath(pkgdir(MIToS), "scripts")

Pkg.activate(SCRIPTS_PATH, io=devnull)
Pkg.instantiate(io=devnull)

using Distributed

using Dates
using MIToS.Utils.Scripts

# ArgParse
# --------

using ArgParse

"""
Parse MIToS scripts command line arguments.
"""
function parse_commandline(
    args...;
    description::AbstractString = "Made with MIToS",
    output::AbstractString = ".mitos.",
    stdout::Bool = true,
    mitos_version = "",
)
    settings = ArgParseSettings(
        description = description,
        version = "MIToS $mitos_version",
        add_version = true,
        epilog = """
                 \n
                 MIToS $mitos_version\n
                 \n
                 Bioinformatics Unit\n
                 Leloir Institute Foundation\n
                 Av. Patricias Argentinas 435, CP C1405BWE, Buenos Aires, Argentina
                 """,
    )

    add_arg_table!(
        settings,
        "FILE",
        Dict(
            :help => "File name. If it is not used, the script reads from STDIN.",
            :required => false,
        ),
        ["--list", "-l"],
        Dict(
            :help => "The input is a list of file names. If -p is used, files will be processed in parallel.",
            :action => :store_true,
        ),
        ["--output", "-o"],
        Dict(
            :help => string(
                """Name of the output file. Output will be gzip if the extension is ".gz".
If it starts with a dot, the name is used as a suffix or extension of the input filename.
If it ends with a dot, is used as a prefix.
If the output name starts and ends with dots, it's used as an interfix before the extension.""",
                stdout ?
                """If a single file is used and there is not a file name (STDIN), the output will be print into
       STDOUT, unless a output filename is used. You can use "STDOUT" to force print into STDOUT.
       STDOUT can not be use with --list.""" : "",
            ),
            :arg_type => AbstractString,
            :default => output,
        ),
        ["--parallel", "-p"],
        Dict(:help => "Number of worker processes.", :arg_type => Int, :default => 1),
        args...,
    )

    return parse_args(settings)
end

# Run script
# ----------

"""
If FILE is not used, calls run_single_script with STDIN. Otherwise calls run_single_script with the FILE.
If `list` is `true`, for loop or `pmap` is used over eachline of the input.
"""
function runscript(args)
    file = args["FILE"]
    if args["list"]
        fh_in = file === nothing ? stdin : open(file, "r")
        if nprocs() != 1
            pmap(eachline(fh_in)) do line
                run_single_script(line, args)
            end
        else
            for line in eachline(fh_in)
                run_single_script(line, args)
            end
        end
    else
        run_single_script(file === nothing ? stdin : file, args)
    end
end

# Distrubuted
# -----------

"""
Adds the needed number of workers.
"""
function set_parallel(parallel)
    N = nprocs()
    if N < parallel
        addprocs(parallel - N + 1)
    end
    nothing
end
