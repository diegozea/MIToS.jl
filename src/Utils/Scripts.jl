"""
Return the version of the loaded module.

Source: https://stackoverflow.com/questions/60587041/julia-getting-the-version-number-of-my-module
"""
loadedversion(m::Module) = VersionNumber(
    Pkg.TOML.parsefile(
        abspath(string(first(methods(m.eval)).file), "..", "..", "Project.toml"),
    )["version"],
)

function _generate_output_name(file, output)
    begins = startswith(output, '.')
    ends = endswith(output, '.')

    if file !== nothing && file !== stdin

        if begins && ends
            parts = split(file, '.')
            return (string(join(parts[1:end-1], '.'), output, parts[end]))
        elseif begins
            return (string(file, output))
        elseif ends
            return (string(output, file))
        else
            return (output)
        end

    elseif !begins && !ends

        return (output)

    else

        return ("STDOUT")

    end
end

"""
Opens the output file or returns STDOUT.
"""
function open_output(file, output)
    output_name = _generate_output_name(file, output)
    if output_name != "STDOUT"
        fh = open(output_name, "w")
        if endswith(output_name, ".gz")
            return (GzipCompressorStream(fh))
        end
        return (fh)
    else
        return (stdout)
    end
end

"""
Close output (check if output is STDOUT).
"""
function close_output(fh_out)
    fh_out === stdout ? nothing : close(fh_out)
    nothing
end

"""
Adds the needed number of workers.
"""
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
Decides if read_file or parse_file, uses parse_file with STDIN
"""
readorparse(input::AbstractString, args...; kargs...) = read_file(input, args...; kargs...)
readorparse(input::Base.LibuvStream, args...; kargs...) =
    parse_file(input, args...; kargs...)
