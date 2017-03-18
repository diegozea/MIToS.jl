const templ = joinpath(splitdir(dirname(@__FILE__))[1], "scripts", "_template.jl")

const list_file = joinpath(dirname(@__FILE__), "data", "list.txt")
const file_one  = joinpath(dirname(@__FILE__), "data", "one.txt")
const file_two  = joinpath(dirname(@__FILE__), "data", "two.txt")

print("""

Test scripts
============
""")

print("""

STDIN -> STDOUT
---------------
""")

let out = readall( pipeline(`$JULIA_HOME/julia -e 'print(trues(2,2))'`, `$JULIA_HOME/julia $templ`) )

    @test ismatch(r"RUN : 0", out)
    @test length(matchall(r"true", out)) == 4
end

@unix_only begin

    print("""
    --arg
    """)

    let out = readall( pipeline(`cat $list_file`, `$JULIA_HOME/julia $templ --arg 42`) )

        @test ismatch(r"RUN : 42", out)
        @test length(matchall(r"\.txt", out)) == 2
    end

    print("""
    --list & -a
    """)

    let out = readall( pipeline(`cat $list_file`, `$JULIA_HOME/julia $templ --a 42 --list`) )

        @test ismatch(r"ONE", out) # Printed into STDOUT
        @test ismatch(r"TWO", out) # Printed into STDOUT
        @test !ismatch(r"RUN : 42", out) # Printed into the file
        @test filesize("one.mitos.txt") != 0 # Created file
        @test filesize("two.mitos.txt") != 0 # Created file
        @test ismatch(r"RUN : 42", readall("one.mitos.txt")) # Printed into the file
        @test ismatch(r"RUN : 42", readall("two.mitos.txt")) # Printed into the file

        # Clean up

        if isfile("one.mitos.txt")
            rm("one.mitos.txt")
        end

        if isfile("two.mitos.txt")
            rm("two.mitos.txt")
        end

    end

end

print("""

File -> File
------------
""")

print("""
    --list & -a
    """)

let out = readall( `$JULIA_HOME/julia $templ $list_file --a 42 --list` )

    @test ismatch(r"ONE", out) # Printed into STDOUT
    @test ismatch(r"TWO", out) # Printed into STDOUT
    @test !ismatch(r"RUN : 42", out) # Printed into the file
    @test filesize("one.mitos.txt") != 0 # Created file
    @test filesize("two.mitos.txt") != 0 # Created file
    @test ismatch(r"RUN : 42", readall("one.mitos.txt")) # Printed into the file
    @test ismatch(r"RUN : 42", readall("two.mitos.txt")) # Printed into the file

    # Clean up

    if isfile("one.mitos.txt")
        rm("one.mitos.txt")
    end

    if isfile("two.mitos.txt")
        rm("two.mitos.txt")
    end

end

print("""

-p (parallel)
-------------
""")

print("""
File -> File (--list)
""")

let out = readall( `$JULIA_HOME/julia $templ $list_file -p 2 --a 42 --list` )

    @test ismatch(r"ONE", out) # Printed into STDOUT
    @test ismatch(r"TWO", out) # Printed into STDOUT
    @test !ismatch(r"RUN : 42", out) # Printed into the file
    @test filesize("one.mitos.txt") != 0 # Created file
    @test filesize("two.mitos.txt") != 0 # Created file
    @test ismatch(r"RUN : 42", readall("one.mitos.txt")) # Printed into the file
    @test ismatch(r"RUN : 42", readall("two.mitos.txt")) # Printed into the file

    # Clean up

    if isfile("one.mitos.txt")
        rm("one.mitos.txt")
    end

    if isfile("two.mitos.txt")
        rm("two.mitos.txt")
    end

end

@unix_only begin

    print("""
    STDIN -> File (--list)
    """)

    let out = readall( pipeline(`cat $list_file`, `$JULIA_HOME/julia $templ -p 2 --a 42 --list` ) )

        @test ismatch(r"ONE", out) # Printed into STDOUT
        @test ismatch(r"TWO", out) # Printed into STDOUT
        @test !ismatch(r"RUN : 42", out) # Printed into the file
        @test filesize("one.mitos.txt") != 0 # Created file
        @test filesize("two.mitos.txt") != 0 # Created file
        @test ismatch(r"RUN : 42", readall("one.mitos.txt")) # Printed into the file
        @test ismatch(r"RUN : 42", readall("two.mitos.txt")) # Printed into the file

        # Clean up

        if isfile("one.mitos.txt")
            rm("one.mitos.txt")
        end

        if isfile("two.mitos.txt")
            rm("two.mitos.txt")
        end

    end

end

print("""

Distances.jl
============
""")

let path_script = joinpath(splitdir(dirname(@__FILE__))[1], "scripts", "Distances.jl"),
    path_file   = joinpath(dirname(@__FILE__), "data", "small.pdb"),
    out_intra  = readall(`$JULIA_HOME/julia $path_script $path_file -o STDOUT`),
    out_inter  = readall(`$JULIA_HOME/julia $path_script $path_file --inter -o STDOUT`)

    @test sum(diff(Float64[ parse(Float64,num) for num in  matchall(r"(\d\.\d+)\n",out_intra) ])) == 0.0
    @test length(matchall(r"\n1,A,", out_intra)) == 1
    @test length(matchall(r"\n1,B,", out_intra)) == 1

    @test sum(Float64[ parse(Float64,num) for num in  matchall(r"(\d\.\d+)\n",out_inter) ]) == 4.0*3.7836265671971385
    @test length(matchall(r"\n1,A,", out_inter)) == 5
    @test length(matchall(r"\n1,B,", out_inter)) == 1

end
