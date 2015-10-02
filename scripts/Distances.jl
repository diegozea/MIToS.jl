#!/usr/bin/env julia

using ArgParse

import MIToS.PDB
@everywhere using MIToS.PDB

function parse_commandline()
    s = ArgParseSettings(description = "Calculates residues distance and writes them into a *.distances.csv file.",
                        version = "MIToS $(Pkg.installed("MIToS"))",
                        add_version = true)

    @add_arg_table s begin
        "--file", "-f"
            help = "Input PDB file"
        "--list", "-l"
            help = "File with a list of input PDB files"
        "--distance", "-d"
            help = "The distance to be calculated, options: All, Heavy, CA, CB"
            arg_type = ASCIIString
            default = "All"
        "--format", "-t"
            help = "Format of the PDB file: It should be pdb or pdbml"
            arg_type = ASCIIString
            default = "pdb"
        "--model", "-m"
            help = "The model to be used, * for all"
            arg_type = ASCIIString
            default = "1"
        "--chain", "-c"
            help = "The chain to be used, * for all"
            arg_type = ASCIIString
            default = "*"
        "--group", "-g"
            help = "Group of atoms to be used, should be ATOM, HETATM or * for all"
            arg_type = ASCIIString
            default = "*"
    end

    s.epilog = """
    \n
    MIToS $(Pkg.installed("MIToS"))\n
    \n
    Bioinformatics Unit\n
    Institute Leloir Foundation\n
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

@everywhere function main(input) # input must be a file
  name, ext = splitext(input)
  fh = open(string(name, ".distances.csv"), "w")
  println(fh, "# MIToS ", Pkg.installed("MIToS"), " Distances.jl ", now())
  println(fh, "# used arguments:")
  for (key, value) in Args
    println(fh, "# \t", key, "\t\t", value)
  end
  println(fh, "model_i,chain_i,group_i,pdbe_i,number_i,name_i,model_j,chain_j,group_j,pdbe_j,number_j,name_j,distance")
  try
    dtyp = ascii(Args["distance"])
    form = ascii(Args["format"])
    if form == "pdb"
      pdb = read(input, PDBFile)
    elseif form == "pdbml"
      pdb = read(input, PDBML)
    else
      throw(ErrorException("--format should be pdb or pdbml."))
    end
    res = residues(pdb, ascii(Args["model"]), ascii(Args["chain"]), ascii(Args["group"]), "*")
    N = length(res)
    for i in 1:(N-1)
      for j in (i+1):N
        @inbounds res1 = res[i]
        @inbounds res2 = res[j]
        dist = distance(res1, res2, criteria=dtyp)
        println(fh, res1.id.model, ",", res1.id.chain, ",", res1.id.group, ",", res1.id.PDBe_number, ",", res1.id.number, ",", res1.id.name, ",",
                res2.id.model, ",", res2.id.chain, ",", res2.id.group, ",", res2.id.PDBe_number, ",", res2.id.number, ",", res2.id.name, ",",
                dist )
      end
    end
  catch err
    println("ERROR: ", input)
    println(err)
  finally
    close(fh)
  end
end

pmap(main, FileList) # Run each file in parallel (with -l)

