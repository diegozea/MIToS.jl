#!/usr/bin/env julia

using ArgParse

import MIToS.MSA
import MIToS.Clustering

@everywhere using MIToS.MSA
@everywhere using MIToS.Clustering

function parse_commandline()
    s = ArgParseSettings(description = """Creates an *.description.csv from a Stockholm file with: the number of columns, sequences, clusters after Hobohm clustering at 62% identity.
    Also the mean, standard deviation and quantiles of: sequence coverage of the MSA, gap percentage and percent identity.""",
                        version = "MIToS $(Pkg.installed("MIToS"))",
                        add_version = true)

    @add_arg_table s begin
        "--file", "-f"
            help = "Input file"
        "--list", "-l"
            help = "File with a list of input files"
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

@everywhere FORMAT = fetch(MIToS.MSA.Stockholm)

@everywhere function _describe(io, input)
  aln = read(input, FORMAT);

  println(io, input, ",", "sequences", ",", "number", ",", "", ",", size(aln, 2))
  println(io, input, ",", "columns",   ",", "number", ",", "", ",", size(aln, 2))

  pid = percentidentity(aln, Float16);
  qid = quantile(pid.list, [0., .25, .5, .75, 1.])

  println(io, input, ",", "percentidentity", ",", "quantile", ",", "0.00", ",", qid[1])
  println(io, input, ",", "percentidentity", ",", "quantile", ",", "0.25", ",", qid[2])
  println(io, input, ",", "percentidentity", ",", "quantile", ",", "0.50", ",", qid[3])
  println(io, input, ",", "percentidentity", ",", "quantile", ",", "0.75", ",", qid[4])

  println(io, input, ",", "percentidentity",   ",", "mean", ",", "", ",", mean(pid))
  println(io, input, ",", "percentidentity",   ",", "std",  ",", "", ",", std(pid))

  cov = coverage(aln);
  qcv = quantile(cov, [0., .25, .5, .75, 1.])

  println(io, input, ",", "coverage", ",", "quantile", ",", "0.00", ",", qcv[1])
  println(io, input, ",", "coverage", ",", "quantile", ",", "0.25", ",", qcv[2])
  println(io, input, ",", "coverage", ",", "quantile", ",", "0.50", ",", qcv[3])
  println(io, input, ",", "coverage", ",", "quantile", ",", "0.75", ",", qcv[4])

  println(io, input, ",", "coverage",   ",", "mean", ",", "", ",", mean(cov))
  println(io, input, ",", "coverage",   ",", "std",  ",", "", ",", std(cov))

  gap = gappercentage(aln, 1);
  qgp = quantile(gap, [0., .25, .5, .75, 1.])

  println(io, input, ",", "gappercentage", ",", "quantile", ",", "0.00", ",", qgp[1])
  println(io, input, ",", "gappercentage", ",", "quantile", ",", "0.25", ",", qgp[2])
  println(io, input, ",", "gappercentage", ",", "quantile", ",", "0.50", ",", qgp[3])
  println(io, input, ",", "gappercentage", ",", "quantile", ",", "0.75", ",", qgp[4])

  println(io, input, ",", "gappercentage",   ",", "mean", ",", "", ",", mean(gap))
  println(io, input, ",", "gappercentage",   ",", "std",  ",", "", ",", std(gap))

  hob = hobohmI(aln, 0.62);
  ncu = getnclusters(hob)

  println(io, input, ",", "clusters",   ",", "number", ",", "", ",", ncu)
end

@everywhere function main(input)
  name, ext = splitext(input)
  fh = open(string(name, ".description.csv"), "w")
  try
    _describe(fh, input)
  catch err
    println("ERROR: ", input)
    println(err)
  finally
    close(fh)
  end
end

pmap(main, FileList) # Run each file in parallel (with -l)
