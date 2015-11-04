#!/usr/bin/env julia

using ArgParse
# TO DO -----------------------------------------------------------------------
import MIToS
@everywhere using MIToS
# -----------------------------------------------------------------------------

function parse_commandline()
# TO DO -----------------------------------------------------------------------
    s = ArgParseSettings(description = "MIToS",
# -----------------------------------------------------------------------------
                        version = "MIToS $(Pkg.installed("MIToS"))",
                        add_version = true)

    @add_arg_table s begin
        "--file", "-f"
            help = "Input file"
        "--list", "-l"
            help = "File with a list of input files"
# TO DO -----------------------------------------------------------------------
        "--arg", "-a"
            help = "Argument"
            arg_type = Int
            default = 0
          # eval_arg = true for arg_type = Bool
# -----------------------------------------------------------------------------
    end

    s.epilog = """
    \n
    MIToS $(Pkg.installed("MIToS"))\n
    \n
    Bioinformatics Unit\n
    Leloir Institute Foundation\n
    Av. Patricias Argentinas 435, CP C1405BWE, Buenos Aires, Argentina
    """

    return parse_args(s)
end

const parsed = parse_commandline()

function _file_names(args)
  file = args["file"]
  list = args["list"]
  if file !== nothing && list === nothing
    return ASCIIString[ file ]
  elseif list !== nothing && file === nothing
    return ASCIIString[ chomp(line) for line in open(readlines, list, "r") ]
  else
    throw(ErrorException("You must use --file or --list and the filename; --help for more information."))
  end
end

const files = _file_names(parsed)

@everywhere Args = remotecall_fetch(1,()->parsed) # Parsed ARGS for each worker
@everywhere FileList = remotecall_fetch(1,()->files) # List of Files for each worker

# TO DO -----------------------------------------------------------------------
@everywhere function main(input) # input must be a file
  try
    arg_one = Args["arg"]
    println("RUN : $arg_one : $input")
  catch err
    println("ERROR: ", input)
    println(err)
  end
end
# -----------------------------------------------------------------------------

pmap(main, FileList) # Run each file in parallel (with -l)
