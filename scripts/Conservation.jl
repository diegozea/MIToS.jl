#!/usr/bin/env julia

using Pkg
using Dates
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
    This takes a MSA file as input and it calculates and saves on *.conservation.csv the
    Shannon entropy (H) and Kullback-Leibler divergence (KL) values for each column
    (Johansson and Toh 2010).\n
    It is possible to do a sequence clustering using the Hobohm I algorithm to avoid the
    effect of sequence redundancy in the conservation scores. Each sequence in a cluster
    is weighted using the inverse of the number of elements in that cluster.\n
    Shannon entropy is a common measure of the residue variability of a particular MSA column.
    For each column, we consider the frequency of the 20 natural protein residues. This uses
    the Euler's number e as the base of the logarithm, so the entropy is measured in nats.\n
    The Kullback-Leibler divergence, also called relative entropy, is a measure of residue
    conservation. It measures how much a probability distribution differs from a background
    distribution. In particular, this implementation measures the divergence between the
    residue distribution of an MSA column and the probabilities derived from the BLOSUM62
    substitution matrix.\n
    Johansson, F., Toh, H., 2010.
    A comparative study of conservation and variation scores.
    BMC Bioinformatics 11, 388.
    """,
    output=".conservation.csv"
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
        println(fh_out, "# MIToS ", Pkg.installed()["MIToS"], " Conservation.jl ", now())
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
        println(fh_out, "# Number of sequences: ", nsequences(msa))
        if args["clustering"]
            clusters = hobohmI(msa, args["threshold"])
            println(fh_out, "# Number of sequence clusters: ", nclusters(clusters))
        else
            clusters = NoClustering()
        end
        probability_table = Probabilities(ContingencyTable(Float64, Val{1}, UngappedAlphabet()))

        KL = mapcolfreq!(kullback_leibler, msa, probability_table, weights = clusters)
        H  = mapcolfreq!(entropy, msa, probability_table, weights = clusters)

        col_names = names(KL,2)
        println(fh_out, "i,H,KL")
        for i in 1:length(col_names)
            println(fh_out, col_names[i], ",", H[i], ",", KL[i])
        end
        # ------------------------------------------------------------------------
    end

end

runscript(args)
