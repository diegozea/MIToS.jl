#!/usr/bin/env julia

using MIToS
include(joinpath(pkgdir(MIToS), "scripts", "_setup_script.jl"))

using CodecZlib
using MIToS.Utils # get_n_words, check_file
using ProgressMeter

function parse_commandline()
    mitos_version = loadedversion(MIToS)
    s = ArgParseSettings(
        description = "Splits a file with multiple sequence alignments in Stockholm format, creating one compressed file per MSA in Stockholm format: accessionumber.gz",
        version = "MIToS $(mitos_version)",
        add_version = true,
    )

    @add_arg_table! s begin
        """
        file
        """
        help = "Input file"
        required = true
        "--path", "-p"
        help = "Path for the output files [default: execution directory]"
        arg_type = String
        default = ""
        """
        --hideprogress
        """
        help = "Hide the progress bar"
        action = :store_true
    end

    s.epilog = """

    \n
 MIToS $(mitos_version)\n
    \n
    Bioinformatics Unit\n
    Leloir Institute Foundation\n
    Av. Patricias Argentinas 435, CP C1405BWE, Buenos Aires, Argentina
    """

    return parse_args(s)
end

const Args = parse_commandline()

function main(input)
    infh = GzipDecompressorStream(open(input))
    lines = []
    id = "no_accessionumber"

    if !Args["hideprogress"]
        totalsize = filesize(input)
        val = 0
        prog = Progress(totalsize, 1)
    end

    # for line in eachline(infh)
    while !eof(infh)
        line = readline(infh)
        if length(line) > 7 && line[1:7] == "#=GF AC"
            id = get_n_words(String(chomp(line)), 3)[3]
        end
        push!(lines, line)
        if startswith(line, "//")
            filename = abspath(joinpath(Args["path"], string(id, ".gz")))
            outfh = GzipCompressorStream(open(filename, "w"))
            for l in lines
                println(outfh, chomp(string(l)))
            end
            close(outfh)
            id = "no_accessionumber"
            empty!(lines)
            if !Args["hideprogress"]
                val += filesize(filename)
                ProgressMeter.update!(prog, val)
            end
        end
    end
    close(infh)

    if !Args["hideprogress"]
        println()
    end
end

main(check_file(abspath(Args["file"])))
