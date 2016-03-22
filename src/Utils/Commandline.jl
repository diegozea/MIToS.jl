"Parse MIToS scripts command line arguments."
function parse_commandline(args...; description::AbstractString="Made with MIToS",
                                    output::AbstractString=".mitos.")
    mitos_version  = Pkg.installed("MIToS")
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
                      :help => """Name of the output file.
                      If it starts with a dot, the name is used as a suffix or extension of the input filename.
                      If it ends with a dot, is used as a suffix.
                      If the output name starts and ends with dots, it's used as an interfix before the extension.
                      If a single file is used and there is not a file name (STDIN), the output will be print into STDOUT.""",
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

    if file !== nothing && file !== STDIN

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

        return("")

    end
end

"Opens the output file or returns STDOUT."
function open_output(file, output)
    output_name = _generate_output_name(file, output)
    if output_name != ""
        fh = open(output_name, "w")
        return(fh)
    else
        return(STDOUT)
    end
end

"Close output (check if output is STDOUT)."
function close_output(fh_out)
    fh_out === STDOUT ? nothing : close(fh_out)
    nothing
end

"""
If FILE is not used, calls main with STDIN. Otherwise calls main with the FILE.
If `list` is `true`, `pmap`/`map` is used over eachline of the input.
"""
function runfun(fun, file::Void, list::Bool)
    if list
        if nprocs() == 1
            map(x->fun(chomp(x)),  eachline(STDIN))
        else
            pmap(x->fun(chomp(x)), eachline(STDIN))
        end
    else
        fun(STDIN)
    end
end

function runfun(fun, file::AbstractString, list::Bool)
    if list
        open(file, "r") do fh
            if nprocs() == 1
                map(x->fun(chomp(x)),  eachline(fh))
            else
                pmap(x->fun(chomp(x)), eachline(fh))
            end
        end
    else
        fun(file)
    end
end

"Adds the needed number of workers."
function set_parallel(parallel)
    N = nprocs()
    if N < parallel
        addprocs(parallel - N + 1)
    end
    nothing
end

