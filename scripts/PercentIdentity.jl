#!/usr/bin/env julia

using MIToS.Utils.Scripts

Args = parse_commandline(
    # TO DO ----------------------------------------------------------------------
    ["--format", "-f"],
    Dict(
        :help => "Format of the MSA: Stockholm, Raw or FASTA",
        :arg_type => ASCIIString,
        :default => "Stockholm"
    ),
    ["--savelist", "-s"],
    Dict(
        :help => "Create and pidlist.csv file with the percentage identity for each pairwise comparison.",
        :action => :store_true
    ),
    # Keywords...
    description="""
    Calculates the percentage identity between all the sequences of an MSA and creates an *.pidstats.csv file with:
    The number of columns and sequences. The mean, standard deviation, median, minimum and maximum values and first and third quantiles of the percentage identity.
    It could also create and pidlist.csv file with the percentage identity for each pairwise comparison.
    """,
    output=".pidstats.csv"
    # ----------------------------------------------------------------------------
    )

set_parallel(Args["parallel"])

@everywhere begin

    const args = remotecall_fetch(1,()->Args)

    import MIToS.Utils.Scripts: script

    # TO DO ----------------------------------------------------------------------
    using MIToS.MSA
    using PairwiseListMatrices
    # ----------------------------------------------------------------------------

    function script(input::Union{Base.LibuvStream,  AbstractString},
                    args,
                    fh_out::Union{Base.LibuvStream, IO})
        # TO DO ------------------------------------------------------------------
        println(fh_out, "# MIToS ", Pkg.installed("MIToS"), " PercentIdentity.jl ", now())
        println(fh_out, "# used arguments:")
        for (key, value) in args
            println(fh_out, "# \t", key, "\t\t", value)
        end
        form = ascii(args["format"])
        if form == "Stockholm"
            msa = read(input, Stockholm)
        elseif form == "FASTA"
            msa = read(input, FASTA)
        elseif form == "Raw"
            msa = read(input, Raw)
        else
            throw(ErrorException("--format should be Stockholm, Raw or FASTA."))
        end
        savelist = args["savelist"]
        println(fh_out, "ncol,nseq,mean,std,min,firstq,median,thirdq,max")
        plm = percentidentity(msa, Float16)
        mean_pid = mean(plm.list)
        min_pid, max_pid = extrema(plm.list)
        first, med, third = quantile(plm.list, Float64[0.25,0.5,0.75])
        println(fh_out, size(msa, 2), ",", size(msa, 1), ",", mean_pid, ",", stdm(plm.list, mean_pid), ",", min_pid, ",", first, ",", med, ",", third, ",", max_pid)
        flush(fh_out)
        if savelist
            writecsv("pidlist.csv", to_table(plm, false))
        end
        # ------------------------------------------------------------------------
    end

end

runscript(args)

