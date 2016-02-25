#!/usr/bin/env julia

using ArgParse

import MIToS.MSA
import MIToS.Information
import PairwiseListMatrices
@everywhere using MIToS.MSA
@everywhere using MIToS.Information
@everywhere using PairwiseListMatrices

function parse_commandline()
    s = ArgParseSettings(description = """Calculates and saves on *.BLMI.csv a Z score and a corrected MI/MIp using BLOSUM62 based pseudo frequencies and clustering of sequences using Hobohm I.""",
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
        "--beta", "-b"
            help = "β for BLOSUM62 pseudo frequencies"
            arg_type = Float64
            default = 8.512
        "--threshold", "-i"
            help = "Percent identity threshold for sequence clustering (Hobohm I)"
            arg_type = Float64
            default = 0.62
        "--maxgap", "-g"
            help = "Maximum fraction of gaps in positions included in calculation"
            arg_type = Float64
            default = 0.5
        "--apc", "-p"
            help = "Use APC correction (MIp)"
            arg_type = Bool
            default = true
            eval_arg = true
        "--samples", "-s"
            help = "Number of samples for Z-score"
            arg_type = Int
            default = 50
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

# TO DO -----------------------------------------------------------------------
@everywhere function main(input) # input must be a file
  name, ext = splitext(input)
  fh = open(string(name, ".BLMI.csv"), "w")
  println(fh, "# MIToS ", Pkg.installed("MIToS"), " BLMI.jl ", now())
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
    zscore, mip = BLMI(msa, beta=Args["beta"], threshold=Args["threshold"],
                       maxgap=Args["maxgap"], apc=Args["apc"], samples=Args["samples"], fixedgaps=Args["fixedgaps"])
    println(fh, "i,j,", Args["apc"] ? "ZBLMIp" : "ZBLMI", ",", Args["apc"] ? "BLMIp" : "BLMI")
    table = hcat(to_table(zscore, false), to_table(mip, false)[:,3])
    writecsv(fh, table)
  catch err
    println("ERROR: ", input)
    println(err)
  finally
    close(fh)
  end
end
# -----------------------------------------------------------------------------

pmap(main, FileList) # Run each file in parallel (with -l)
