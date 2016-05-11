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
    ["--beta", "-b"],
    Dict(
        :help => "β for BLOSUM62 pseudo frequencies",
        :arg_type => Float64,
        :default => 8.512
    ),
    ["--threshold", "-i"],
    Dict(
        :help => "Percent identity threshold for sequence clustering (Hobohm I)",
        :arg_type => Float64,
        :default => 62.0
    ),
    ["--maxgap", "-g"],
    Dict(
        :help => "Maximum fraction of gaps in positions included in calculation",
        :arg_type => Float64,
        :default => 0.5
    ),
    ["--apc", "-a"],
    Dict(
        :help => "Use APC correction (MIp)",
        :action => :store_false
    ),
    ["--samples", "-s"],
    Dict(
        :help => "Number of samples for Z-score",
        :arg_type => Int,
        :default => 50
    ),
    ["--fixedgaps", "-F"],
    Dict(
        :help => "Fix gaps positions for the random samples",
        :action => :store_false
    ),
    # Keywords...
    description="""
    This takes a MSA file as input.
    Calculates and saves on *.BLMI.csv a Z score and a corrected MI/MIp.
    The script uses BLOSUM62 based pseudo frequencies and sequences clustering (Hobohm I).
    """,
    output=".BLMI.csv"
    # ----------------------------------------------------------------------------
    )

set_parallel(Args["parallel"])

@everywhere begin

    const args = remotecall_fetch(1,()->Args)

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
        println(fh_out, "# MIToS ", Pkg.installed("MIToS"), " BLMI.jl ", now())
        println(fh_out, "# used arguments:")
        for (key, value) in args
            println(fh_out, "# \t", key, "\t\t", value)
        end
        form = ascii(args["format"])
        if form == "Stockholm"
            msa = readorparse(input, Stockholm)
        elseif form == "FASTA"
            msa = readorparse(input, FASTA)
        elseif form == "Raw"
            msa = readorparse(input, Raw)
        else
            throw(ErrorException("--format should be Stockholm, Raw or FASTA."))
        end
        zscore, mip = BLMI(msa, beta=args["beta"], threshold=args["threshold"],
                           maxgap=args["maxgap"], apc=args["apc"], samples=args["samples"], fixedgaps=args["fixedgaps"])
        println(fh_out, "i,j,", args["apc"] ? "ZBLMIp" : "ZBLMI", ",", args["apc"] ? "BLMIp" : "BLMI")
        table = hcat(to_table(zscore, false), to_table(mip, false)[:,3])
        writecsv(fh_out, table)
        # ------------------------------------------------------------------------
    end

end

runscript(args)
