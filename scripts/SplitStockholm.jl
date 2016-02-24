#!/usr/bin/env julia

using ArgParse
using GZip
using MIToS.Utils # get_n_words, check_file

function parse_commandline()
    s = ArgParseSettings(description = "Splits a file with multiple sequence alignments in Stockholm format, creating one compressed file per MSA in Stockholm format: accessionumber.gz",
                         version = "MIToS $(Pkg.installed("MIToS"))",
                         add_version = true)

    @add_arg_table s begin
        "file"
            help = "Input file"
            required = true
        "--path", "-p"
        		help = "Path for the output files [default: execution directory]"
		        arg_type = ASCIIString
    		    default = ""
    end

    s.epilog = """

    \n
    MIToS $(Pkg.installed("MIToS"))\n
    \n
    Bioinformatics Unit\n
    Leloir Institute Foundation\n
    Av. Patricias Argentinas 435, CP C1405BWE, Buenos Aires, Argentina
    """

    return parse_args(s)
end

const Args = parse_commandline()

function main(input)
	infh = GZip.open(input)
	lines = []
	id = "no_accessionumber"
	for line in eachline(infh)
		if length(line) > 7 && line[1:7] == "#=GF AC"
			id = get_n_words(line, 3)[3]
		end
		push!(lines, line)
		if line == "//\n"
			filename = joinpath(Args["path"], string(id, ".gz"))
			outfh = GZip.open(filename, "w")
			for l in lines
				write(outfh, string(l))
			end
			close(outfh)
			id = "no_accessionumber"
			empty!(lines)
		end
	end
	close(infh)
end

main(check_file(Args["file"]))
