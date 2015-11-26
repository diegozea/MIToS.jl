#!/usr/bin/env julia

using ArgParse

import MIToS.MSA
@everywhere using MIToS.MSA
import MIToS.Clustering
@everywhere using MIToS.Clustering
import PairwiseListMatrices
@everywhere using PairwiseListMatrices

function parse_commandline()
    s = ArgParseSettings(description = """Calculates the percentage identity between all the sequences of an MSA and creates an *.pidstats.csv file with:
    mean, standard deviation, median, minimum and maximum values of the percentage identity and the number of columns and sequences.
    It could also create and *.pidlist.csv file with the percentage identity for each pairwise comparison.""",
                        version = "MIToS $(Pkg.installed("MIToS"))",
                        add_version = true)

    @add_arg_table s begin
        "--file", "-f"
            help = "Input file"
        "--list", "-l"
            help = "File with a list of input files"
        "--format", "-t"
            help = "Format of the MSA: stockholm, raw or fasta"
            arg_type = ASCIIString
            default = "stockholm"
        "--pidlist", "-p"
            help = "Create and *.pidlist.csv file with the percentage identity for each pairwise comparison."
            arg_type = Bool
            default = false
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

@everywhere function main(input)
  name, ext = splitext(input)
  savelist = Args["pidlist"]
  fh = open(string(name, ".pidstats.csv"), "w")
  println(fh, "# MIToS ", Pkg.installed("MIToS"), " PercentIdentity.jl ", now())
  println(fh, "# used arguments:")
  for (key, value) in Args
    println(fh, "# \t", key, "\t\t", value)
  end
  println(fh, "ncol,nseq,mean,std,min,max,median")
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
    plm = percentidentity(msa, Float16)
    println(fh, size(msa, 2), ",", size(msa, 1), ",", mean(plm.list), ",", std(plm.list), ",", minimum(plm.list), ",", maximum(plm.list), ",", median(plm.list))
    flush(fh)
    savelist && writecsv(string(name, ".pidlist.csv"), to_table(plm, false))
  catch err
    println("ERROR: ", input)
    println(err)
  finally
    close(fh)
  end
end

pmap(main, FileList) # Run each file in parallel (with -l)
