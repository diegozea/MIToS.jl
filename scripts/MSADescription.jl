#!/usr/bin/env julia

using MIToS.Utils.Scripts

Args = parse_commandline(
    ["--format", "-f"],
    Dict(
        :help => "Format of the MSA: Stockholm, Raw or FASTA",
        :arg_type => ASCIIString,
        :default => "Stockholm"
    ),
    description="""
    Creates an *.description.csv from a Stockholm file with: the number of columns, sequences, clusters after Hobohm clustering at 62% identity and mean percent identity.
    Also the mean, standard deviation and quantiles of: sequence coverage of the MSA, gap percentage.
    """,
    output=".description.csv"
    # ----------------------------------------------------------------------------
    )

set_parallel(Args["parallel"])

@everywhere begin

    const args = remotecall_fetch(1,()->Args)

    import MIToS.Utils.Scripts: script

    # TO DO ----------------------------------------------------------------------
    using MIToS.MSA

    function _describe(fh_out, input, format)

    end
    # ----------------------------------------------------------------------------

    function script(input::Union{Base.LibuvStream,  AbstractString},
                    args,
                    fh_out::Union{Base.LibuvStream, IO})
        # TO DO ------------------------------------------------------------------
        println(fh_out, "# MIToS ", Pkg.installed("MIToS"), " MSADescription.jl ", now())
        println(fh_out, "# used arguments:")

        for (key, value) in args
            println(fh_out, "# \t", key, "\t\t", value)
        end

        form = ascii(args["format"])

        if form == "Stockholm"
            aln = read(input, Stockholm)
        elseif form == "FASTA"
            aln = read(input, FASTA)
        elseif form == "Raw"
            aln = read(input, Raw)
        else
            throw(ErrorException("--format should be Stockholm, Raw or FASTA."))
        end

        println(fh_out, input, ",", "sequences", ",", "number", ",", "", ",", size(aln, 1))
        println(fh_out, input, ",", "columns",   ",", "number", ",", "", ",", size(aln, 2))

        cov = coverage(aln);
        qcv = quantile(cov, [0., .25, .5, .75, 1.])

        println(fh_out, input, ",", "coverage", ",", "quantile", ",", "0.00", ",", qcv[1])
        println(fh_out, input, ",", "coverage", ",", "quantile", ",", "0.25", ",", qcv[2])
        println(fh_out, input, ",", "coverage", ",", "quantile", ",", "0.50", ",", qcv[3])
        println(fh_out, input, ",", "coverage", ",", "quantile", ",", "0.75", ",", qcv[4])

        println(fh_out, input, ",", "coverage",   ",", "mean", ",", "", ",", mean(cov))
        println(fh_out, input, ",", "coverage",   ",", "std",  ",", "", ",", std(cov))

        println(fh_out, input, ",", "percentidentity",   ",", "mean", ",", "", ",", meanpercentidentity(aln))

        gap = gapfraction(aln, 1);
        qgp = quantile(gap, [0., .25, .5, .75, 1.])

        println(fh_out, input, ",", "gapfraction", ",", "quantile", ",", "0.00", ",", qgp[1])
        println(fh_out, input, ",", "gapfraction", ",", "quantile", ",", "0.25", ",", qgp[2])
        println(fh_out, input, ",", "gapfraction", ",", "quantile", ",", "0.50", ",", qgp[3])
        println(fh_out, input, ",", "gapfraction", ",", "quantile", ",", "0.75", ",", qgp[4])

        println(fh_out, input, ",", "gapfraction",   ",", "mean", ",", "", ",", mean(gap))
        println(fh_out, input, ",", "gapfraction",   ",", "std",  ",", "", ",", std(gap))

        hob = hobohmI(aln, 62);
        ncu = nclusters(hob)

        println(fh_out, input, ",", "clusters",   ",", "number", ",", "", ",", ncu)
        # ------------------------------------------------------------------------
    end

end

runscript(args)
