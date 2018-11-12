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
    ["--lambda", "-L"],
    Dict(
        :help => "Low count value",
        :arg_type => Float64,
        :default => 0.05
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
        :default => 100
    ),
    ["--usegap", "-G"],
    Dict(
        :help => "Use gaps on statistics",
        :action => :store_true
    ),
    ["--fixedgaps", "-F"],
    Dict(
        :help => "Fix gaps positions for the random samples",
        :action => :store_false
    ),
    # Keywords...
    description="""
    This takes a MSA file as input.
    It calculates and saves on *.busjle09.csv a Z score and a corrected MI/MIp as described on:\n
    Buslje, C. M., Santos, J., Delfino, J. M., & Nielsen, M. (2009).
    Correction for phylogeny, small number of observations and data redundancy improves the identification of coevolving amino acid pairs using mutual information.
    Bioinformatics, 25(9), 1125-1131.
    """,
    output=".busjle09.csv"
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
        println(fh_out, "# MIToS ", Pkg.installed()["MIToS"], " Buslje09.jl ", now())
        println(fh_out, "# \tBuslje, C. M., Santos, J., Delfino, J. M., & Nielsen, M. (2009).")
        println(fh_out, "# \tCorrection for phylogeny, small number of observations and data redundancy improves the identification of coevolving amino acid pairs using mutual information.")
        println(fh_out, "# \tBioinformatics, 25(9), 1125-1131.")
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
        zscore, mip = buslje09(msa, lambda=args["lambda"],
                               clustering=args["clustering"], threshold=args["threshold"],
                               maxgap=args["maxgap"], apc=args["apc"], samples=args["samples"],
                               alphabet = args["usegap"] ? GappedAlphabet() : UngappedAlphabet(),
                               fixedgaps=args["fixedgaps"])
        println(fh_out, "i,j,", args["apc"] ? "ZMIp" : "ZMI", ",", args["apc"] ? "MIp" : "MI")
        table = hcat(to_table(zscore, diagonal=false), to_table(mip, diagonal=false)[:,3])
        writedlm(fh_out, table, ',')
        # ------------------------------------------------------------------------
    end

end

runscript(args)
