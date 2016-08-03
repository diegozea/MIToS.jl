#!/usr/bin/env julia

using MIToS.Utils.Scripts

Args = parse_commandline(
    # TO DO ----------------------------------------------------------------------
    ["--distance", "-d"],
    Dict(
        :help => "The distance to be calculated, options: All, Heavy, CA, CB",
        :arg_type => ASCIIString,
        :default => "All"
    ),
    ["--format", "-f"],
    Dict(
        :help => "Format of the PDB file: It should be PDBFile or PDBML",
        :arg_type => ASCIIString,
        :default => "PDBFile"
    ),
    ["--model", "-m"],
    Dict(
        :help => "The model to be used, * for all",
        :arg_type => ASCIIString,
        :default => "1"
    ),
    ["--chain", "-c"],
    Dict(
        :help => "The chain to be used, * for all",
        :arg_type => ASCIIString,
        :default => "*"
    ),
    ["--group", "-g"],
    Dict(
        :help => "Group of atoms to be used, should be ATOM, HETATM or * for all",
        :arg_type => ASCIIString,
        :default => "*"
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

    const args = remotecall_fetch(1,()->Args)

    import MIToS.Utils.Scripts: script

    # TO DO ----------------------------------------------------------------------
    using MIToS.PDB
    # ----------------------------------------------------------------------------

    function script(input::Union{Base.LibuvStream,  AbstractString},
                    args,
                    fh_out::Union{Base.LibuvStream, IO})
        # TO DO ------------------------------------------------------------------
        println(fh_out, "# MIToS ", Pkg.installed("MIToS"), " Distances.jl ", now())
        println(fh_out, "# used arguments:")
        for (key, value) in args
            println(fh_out, "# \t", key, "\t\t", value)
        end
        println(fh_out, "model_i,chain_i,group_i,pdbe_i,number_i,name_i,model_j,chain_j,group_j,pdbe_j,number_j,name_j,distance")
        dtyp = ascii(args["distance"])
        form = ascii(args["format"])
        if form == "PDBFile"
            res = readorparse(input, PDBFile)
        elseif form == "PDBML"
            res = readorparse(input, PDBML)
        else
            throw(ErrorException("--format should be PDBFile or PDBML."))
        end
        res = residues(res, ascii(args["model"]), ascii(args["chain"]), ascii(args["group"]), "*")
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
