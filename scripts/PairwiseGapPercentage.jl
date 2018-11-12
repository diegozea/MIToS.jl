#!/usr/bin/env julia

using Pkg
using Dates
using DelimitedFiles
using Distributed
using MIToS.Utils.Scripts

Args = parse_commandline(
    # TO DO ----------------------------------------------------------------------
    ["--format", "-f"],
    Dict(
        :help => "Format of the MSA: Stockholm, Raw or FASTA",
        :arg_type => String,
        :default => "Stockholm"
    ),
    ["--clustering", "-c"],
    Dict(
        :help => "Sequence clustering (Hobohm I)",
        :action => :store_false
    ),
    ["--threshold", "-i"],
    Dict(
        :help => "Percent identity threshold for sequence clustering (Hobohm I)",
        :arg_type => Float64,
        :default => 62.0
    ),
    # Keywords...
    description="""
    Calculates and saves on *.pairwisegaps.csv the percentage of gaps on columns pairs (union and intersection) using sequence clustering (Hobohm I).
    """,
    output=".pairwisegaps.csv"
    # ----------------------------------------------------------------------------
    )

set_parallel(Args["parallel"])

@everywhere begin

    const args = remotecall_fetch(()->Args,1)

    import MIToS.Utils.Scripts: script

    # TO DO ----------------------------------------------------------------------
    using MIToS.MSA
    using MIToS.Information
    using PairwiseListMatrices
    # ----------------------------------------------------------------------------

    function script(input::Union{Base.LibuvStream,  AbstractString},
                    args,
                    fh_out::Union{Base.LibuvStream, IO})
        # TO DO ------------------------------------------------------------------
        println(fh_out, "# MIToS ", Pkg.installed()["MIToS"], " PairwiseGapPercentage.jl ", now())
        println(fh_out, "# used arguments:")
        for (key, value) in args
            println(fh_out, "# \t", key, "\t\t", value)
        end
        form = string(args["format"])
        if form == "Stockholm"
            msa = readorparse(input, Stockholm)
        elseif form == "FASTA"
            msa = readorparse(input, FASTA)
        elseif form == "Raw"
            msa = readorparse(input, Raw)
        else
            throw(ErrorException("--format should be Stockholm, Raw or FASTA."))
        end
        gapsunion, gapsinter = pairwisegapfraction(msa, clustering=args["clustering"], threshold=args["threshold"])
        println(fh_out, "i,j,gapunion,gapintersection")
        table = hcat(to_table(gapsunion), to_table(gapsinter)[:,3])
        writedlm(fh_out, table, ',')
        # ------------------------------------------------------------------------
    end

end

runscript(args)
