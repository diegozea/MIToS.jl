#!/usr/bin/env julia

using ArgParse

import MIToS.MSA
import MIToS.Information
import PairwiseListMatrices
@everywhere using MIToS.MSA
@everywhere using MIToS.Information
@everywhere using PairwiseListMatrices

function parse_commandline()
    s = ArgParseSettings(description = """Calculates and saves on *.busjle09.csv a Z score and a corrected MI/MIp as described on:\n\n\n
    Buslje, C. M., Santos, J., Delfino, J. M., & Nielsen, M. (2009). Correction for phylogeny, small number of observations and data redundancy improves the identification of coevolving amino acid pairs using mutual information. Bioinformatics, 25(9), 1125-1131.""",
                        version = "MIToS $(Pkg.installed("MIToS"))",
                        add_version = true)

    @add_arg_table s begin
        "--file", "-f"
            help = "Input MSA file"
        "--list", "-l"
            help = "File with a list of input MSA files"
        "--format", "-t"
            help = "Format of the MSA: stockholm, raw or fasta"
            arg_type = ASCIIString
            default = "stockholm"
        "--lambda", "-a"
            help = "Low count value"
            arg_type = Float64
            default = 0.05
        "--clustering", "-c"
            help = "Sequence clustering (Hobohm I)"
            arg_type = Bool
            default = true
            eval_arg = true
        "--threshold", "-i"
            help = "Percent identity threshold for clustering"
            arg_type = Float64
            default = 62.0
        "--maxgap", "-g"
            help = "Maximum fraction of gaps in positions included in calculation"
            arg_type = Float64
            default = 0.5
        "--apc", "-p"
            help = "Use APC correction (MIp)"
            arg_type = Bool
            default = true
            eval_arg = true
        "--usegap", "-G"
            help = "Use gaps on statistics"
            arg_type = Bool
            default = false
            eval_arg = true
        "--samples", "-s"
            help = "Number of samples for Z-score"
            arg_type = Int
            default = 100
        "--fixedgaps", "-F"
            help = "Fix gaps positions for the random samples"
            arg_type = Bool
            default = true
            eval_arg = true
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

@everywhere function main(input) # input must be a file
  name, ext = splitext(input)
  fh = open(string(name, ".buslje09.csv"), "w")
  println(fh, "# MIToS ", Pkg.installed("MIToS"), " Buslje09.jl ", now())
  println(fh, "# \tBuslje, C. M., Santos, J., Delfino, J. M., & Nielsen, M. (2009).")
  println(fh, "# \tCorrection for phylogeny, small number of observations and data redundancy improves the identification of coevolving amino acid pairs using mutual information.")
  println(fh, "# \tBioinformatics, 25(9), 1125-1131.""")
  println(fh, "# used arguments:")
  for (key, value) in Args
    println(fh, "# \t", key, "\t\t", value)
  end
  try
    form = ascii(Args["format"])
    if form == "stockholm"
      msa = read(input, Stockholm)
    elseif form == "fasta"
      msa = read(input, FASTA)
    elseif form == "raw"
      msa = read(input, Raw)
    else
      throw(ErrorException("--format should be stockholm, raw or fasta."))
    end
    zscore, mip = buslje09(msa,lambda=Args["lambda"],
                               clustering=Args["clustering"], threshold=Args["threshold"],
                               maxgap=Args["maxgap"], apc=Args["apc"], samples=Args["samples"],
                               usegap=Args["usegap"], fixedgaps=Args["fixedgaps"])
    println(fh, "i,j,", Args["apc"] ? "ZMIp" : "ZMI", ",", Args["apc"] ? "MIp" : "MI")
    table = hcat(to_table(zscore, false), to_table(mip, false)[:,3])
    writecsv(fh, table)
  catch err
    println("ERROR: ", input)
    println(err)
  finally
    close(fh)
  end
end

pmap(main, FileList) # Run each file in parallel (with -l)
