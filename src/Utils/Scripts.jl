"Parse MIToS scripts command line arguments."
function parse_commandline(args...; description::AbstractString="Made with MIToS",
                                    output::AbstractString=".mitos.",
                                    stdout::Bool=true)
    mitos_version  = Pkg.installed()["MIToS"]
    settings = ArgParseSettings(description = description,
                                version = "MIToS $mitos_version",
                                add_version = true,
                                epilog =    """
                                            \n
                                            MIToS $mitos_version\n
                                            \n
                                            Bioinformatics Unit\n
                                            Leloir Institute Foundation\n
                                            Av. Patricias Argentinas 435, CP C1405BWE, Buenos Aires, Argentina
                                            """
                                )

    add_arg_table(settings,
                  "FILE",
                  Dict(
                      :help => "File name. If it is not used, the script reads from STDIN.",
                      :required => false
                      ),
                  ["--list", "-l"],
                  Dict(
                      :help => "The input is a list of file names. If -p is used, files will be processed in parallel.",
                      :action => :store_true
                      ),
                  ["--output", "-o"],
                  Dict(
                      :help => string( """Name of the output file. Output will be gzip if the extension is ".gz".
                      If it starts with a dot, the name is used as a suffix or extension of the input filename.
                      If it ends with a dot, is used as a prefix.
                      If the output name starts and ends with dots, it's used as an interfix before the extension.""",
                      stdout ? """If a single file is used and there is not a file name (STDIN), the output will be print into
                      STDOUT, unless a output filename is used. You can use "STDOUT" to force print into STDOUT.
                      STDOUT can not be use with --list.""" : ""),
                      :arg_type => AbstractString,
                      :default => output
                      ),
                  ["--parallel", "-p"],
                  Dict(
                      :help => "Number of worker processes.",
                      :arg_type => Int,
                      :default => 1
                      ),

                  args...)

    return parse_args(settings)
end

function _generate_output_name(file, output)
    begins = startswith(output, '.')
    ends = endswith(output,'.')

    if file !== nothing && file !== stdin

        if begins && ends
            parts = split(file, '.')
            return( string(join(parts[1:end-1],'.'), output, parts[end]) )
        elseif begins
            return( string(file, output) )
        elseif ends
            return( string(output, file) )
        else
            return(output)
        end

    elseif !begins && !ends

        return(output)

    else

        return("STDOUT")

    end
end

"Opens the output file or returns STDOUT."
function open_output(file, output)
    output_name = _generate_output_name(file, output)
    if output_name != "STDOUT"
        fh = open(output_name, "w")
        if endswith(output_name, ".gz")
            return(GzipCompressorStream(fh))
        end
        return(fh)
    else
        return(stdout)
    end
end

"Close output (check if output is STDOUT)."
function close_output(fh_out)
    fh_out === stdout ? nothing : close(fh_out)
    nothing
end

"Adds the needed number of workers."
function set_parallel(parallel)
    N = nprocs()
    if N < parallel
        addprocs(parallel - N + 1)
    end
    nothing
end

script(x, y, z) = throw("Define your script function!")

"""
Opens and closes the output stream, tries to run the defined script.
"""
function run_single_script(input, args)
    fh_out = open_output(input, args["output"])
    try
        script(input, args, fh_out)
    catch err
        println("ERROR: ", input)
        println(err)
    finally
        close(fh_out)
    end
end

"""
If FILE is not used, calls run_single_script with STDIN. Otherwise calls run_single_script with the FILE.
If `list` is `true`, for loop or `pmap` is used over eachline of the input.
"""
function runscript(args)
    file = args["FILE"]
    if args["list"]
        fh_in = file === nothing ? stdin : open(file, "r")
        if nprocs() != 1
            pmap(eachline(fh_in)) do line
                run_single_script(line, args)
            end
        else
            for line in eachline(fh_in)
                run_single_script(line, args)
            end
        end
    else
        run_single_script(file === nothing ? stdin : file, args)
    end
end

"Decides if read or parse, uses parse with STDIN"
readorparse(input::AbstractString, args...; kargs...) = read(input, args...; kargs...)
readorparse(input::Base.LibuvStream, args...; kargs...) = parse(input, args...; kargs...)
