#!/usr/bin/env julia

using Pkg
using Dates
using Distributed
using MIToS.Utils.Scripts

Args = parse_commandline(
    # TO DO ----------------------------------------------------------------------
    ["--distance", "-d"],
    Dict(
        :help => "The distance to be calculated, options: All, Heavy, CA, CB",
        :arg_type => String,
        :default => "All"
    ),
    ["--format", "-f"],
    Dict(
        :help => "Format of the PDB file: It should be PDBFile or PDBML",
        :arg_type => String,
        :default => "PDBFile"
    ),
    ["--model", "-m"],
    Dict(
        :help => "The model to be used, use All for all",
        :arg_type => String,
        :default => "1"
    ),
    ["--chain", "-c"],
    Dict(
        :help => "The chain to be used, use All for all",
        :arg_type => String,
        :default => "All"
    ),
    ["--group", "-g"],
    Dict(
        :help => "Group of atoms to be used, should be ATOM, HETATM or All for all",
        :arg_type => String,
        :default => "All"
    ),
    ["--inter", "-i"],
    Dict(
        :help => "Calculate inter chain distances",
        :action => :store_true
    ),
    # Keywords...
    description="""
    Calculates residues distance and writes them into a *.distances.csv.gz gzipped file.
    """,
    output=".distances.csv.gz"
    # ----------------------------------------------------------------------------
    )

set_parallel(Args["parallel"])

@everywhere begin

    const args = remotecall_fetch(()->Args,1)

    import MIToS.Utils.Scripts: script

    # TO DO ----------------------------------------------------------------------
    using MIToS.PDB
    # ----------------------------------------------------------------------------

    function script(input::Union{Base.LibuvStream,  AbstractString},
                    args,
                    fh_out::Union{Base.LibuvStream, IO})
        # TO DO ------------------------------------------------------------------
        println(fh_out, "# MIToS ", Pkg.installed()["MIToS"], " Distances.jl ", now())
        println(fh_out, "# used arguments:")
        for (key, value) in args
            println(fh_out, "# \t", key, "\t\t", value)
        end
        println(fh_out, "model_i,chain_i,group_i,pdbe_i,number_i,name_i,model_j,chain_j,group_j,pdbe_j,number_j,name_j,distance")
        dtyp = string(args["distance"])
        form = string(args["format"])
        if form == "PDBFile"
            res = readorparse(input, PDBFile)
        elseif form == "PDBML"
            res = readorparse(input, PDBML)
        else
            throw(ErrorException("--format should be PDBFile or PDBML."))
        end
        model_arg = string(args["model"]) == "All" ? All : string(args["model"])
        chain_arg = string(args["chain"]) == "All" ? All : string(args["chain"])
        group_arg = string(args["group"]) == "All" ? All : string(args["group"])
        res = residues(res, model_arg, chain_arg, group_arg, All)
        N = length(res)
        inter = !Bool(args["inter"])
        for i in 1:(N-1)
            for j in (i+1):N
                @inbounds res1 = res[i]
                @inbounds res2 = res[j]
                if inter && res1.id.chain != res2.id.chain
                    continue
                end
                dist = distance(res1, res2, criteria=dtyp)
                println(fh_out, res1.id.model, ",", res1.id.chain, ",", res1.id.group, ",", res1.id.PDBe_number, ",", res1.id.number, ",", res1.id.name, ",",
                        res2.id.model, ",", res2.id.chain, ",", res2.id.group, ",", res2.id.PDBe_number, ",", res2.id.number, ",", res2.id.name, ",",
                        dist )
            end
        end
        # ------------------------------------------------------------------------
    end

end

runscript(args)
