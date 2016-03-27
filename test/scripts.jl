const templ = joinpath(Pkg.dir("MIToS"), "scripts", "_template.jl")

const list_file = joinpath(Pkg.dir("MIToS"), "test", "data", "list.txt")
const file_one  = joinpath(Pkg.dir("MIToS"), "test", "data", "one.txt")
const file_two  = joinpath(Pkg.dir("MIToS"), "test", "data", "two.txt")

print("""

Test scripts
============
""")

print("""

STDIN -> STDOUT
---------------
""")

let out = readall( pipeline(`julia -e 'print(trues(2,2))'`, `julia $templ`) )

    @test ismatch(r"RUN : 0", out)
    @test length(matchall(r"true", out)) == 4
end

@unix_only begin

    print("""
    --arg
    """)

    let out = readall( pipeline(`cat $list_file`, `julia $templ --arg 42`) )

        @test ismatch(r"RUN : 42", out)
        @test length(matchall(r"\.txt", out)) == 2
    end

    print("""
    --list & -a
    """)

    let out = readall( pipeline(`cat $list_file`, `julia $templ --a 42 --list`) )

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

let out = readall( `julia $templ $list_file --a 42 --list` )

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

let out = readall( `julia $templ $list_file -p 2 --a 42 --list` )

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

    let out = readall( pipeline(`cat $list_file`, `julia $templ -p 2 --a 42 --list` ) )

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
