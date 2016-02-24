#!/usr/bin/env julia

using ArgParse
using MIToS.MSA

function parse_commandline()
    s = ArgParseSettings(description = """Creates a file in Stockholm format with the aligned columns from a Pfam Stockholm file.
    Insertions are deleted, as they are unaligned in a proﬁle HMM.
    The output file *.aligned.* contains as annotations UniProt residue numbers and column numbers in the original MSA.""",
                        version = "MIToS $(Pkg.installed("MIToS"))",
                        add_version = true)

    @add_arg_table s begin
        "--file", "-f"
            help = "Pfam stockholm file"
        "--list", "-l"
            help = "File with a list of Pfam stockholm files"
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

@everywhere FORMAT = fetch(MIToS.MSA.Stockholm)

@everywhere function main(input)
  try
    name, ext = splitext(input)
    aln = read(input, FORMAT, generatemapping=true, useidcoordinates=true, deletefullgaps=true)
    write(string(name, ".aligned", ext), aln, FORMAT)
  catch err
    println("ERROR: ", input)
    println(err)
  end
end

pmap(main, FileList) # Run each file in parallel (with -l)
