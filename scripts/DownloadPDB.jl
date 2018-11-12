#!/usr/bin/env julia

using Pkg
using Dates
using Distributed
using ArgParse
import MIToS
@everywhere using MIToS

function parse_commandline()
    s = ArgParseSettings(description = "Download gzipped files from PDB.",
                        version = "MIToS $(Pkg.installed()["MIToS"])",
                        add_version = true)

    @add_arg_table s begin
        "--code", "-c"
            help = "PDB code"
        "--list", "-l"
            help = "File with a list of PDB codes (one per line)"
        "--format", "-f"
            help = "Format. It should be PDBFile (pdb) or PDBML (xml)"
            arg_type = String
            default = "PDBML"
    end

    s.epilog = """
    \n
    MIToS $(Pkg.installed()["MIToS"])\n
    \n
    Bioinformatics Unit\n
    Leloir Institute Foundation\n
    Av. Patricias Argentinas 435, CP C1405BWE, Buenos Aires, Argentina
    """

    return parse_args(s)
end

const parsed = parse_commandline()

function _file_names(args)
  file = args["code"]
  list = args["list"]
  if file !== nothing && list === nothing
    return String[ file ]
  elseif list !== nothing && file === nothing
    return String[ chomp(line) for line in open(readlines, list, "r") ]
  else
    throw(ErrorException("You must use --code or --list and the filename; --help for more information."))
  end
end

const files = _file_names(parsed)

@everywhere Args = remotecall_fetch(()->parsed,1) # Parsed ARGS for each worker
@everywhere FileList = remotecall_fetch(()->files,1) # List of Files for each worker

@everywhere function main(input) # input must be a file
  try
    format = Args["format"]
    if format == "PDBML" || format == "PDBFile"
      MIToS.PDB.downloadpdb(input, format=(format == "PDBML" ? MIToS.PDB.PDBML : MIToS.PDB.PDBFile))
    else
      throw(ErrorException("--format should be PDBML or PDBFile"))
    end
  catch err
    println("ERROR: ", input)
    println(err)
  end
end

pmap(main, FileList) # Run each file in parallel (with -l)
