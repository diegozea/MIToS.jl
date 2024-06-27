# Set up function
function _set_up()
    path = tempdir()
    one_file = joinpath(path, "one.tmp")
    two_file = joinpath(path, "two.tmp")
    list_file = joinpath(path, "list.tmp")
    write(one_file, "ONE")
    write(two_file, "TWO")
    write(list_file, one_file, "\n", two_file, "\n")
    one_file, two_file, list_file
end

# Clean up functions
function _delete_file(path)
    if isfile(path)
        rm(path)
    end
end

function _clean_up_outputs()
    path = tempdir()
    _delete_file(joinpath(path, "one.mitos.tmp"))
    _delete_file(joinpath(path, "two.mitos.tmp"))
end

function _clean_up()
    path = tempdir()
    _delete_file(joinpath(path, "list.tmp"))
    _delete_file(joinpath(path, "one.tmp"))
    _delete_file(joinpath(path, "two.tmp"))
end

# Tests
@testset "Template" begin

    # julia bin
    julia = joinpath(Base.Sys.BINDIR, "julia")
    # ../../
    mitos_folder = splitdir(splitdir(dirname(@__FILE__))[1])[1]

    template = joinpath(mitos_folder, "scripts", "_template.jl")
    file_one, file_two, list_file = _set_up()
    one_out = joinpath(tempdir(), "one.mitos.tmp")
    two_out = joinpath(tempdir(), "two.mitos.tmp")

    @testset "STDIN -> STDOUT" begin
        out = read(pipeline(`$julia -e 'print("Hello hello")'`, `$julia $template`), String)

        @test occursin(r"RUN : 0", out)
        @test length(collect((m.match for m in eachmatch(r"[Hh]ello", out)))) == 2

        @static if Sys.isunix()

            @testset "--arg" begin
                out = read(pipeline(`cat $list_file`, `$julia $template --arg 42`), String)

                @test occursin(r"RUN : 42", out)
                @test length(collect((m.match for m in eachmatch(r"\.tmp", out)))) == 2
            end

            @testset "--list & -a" begin
                _clean_up_outputs()
                out = read(
                    pipeline(`cat $list_file`, `$julia $template -a 42 --list`),
                    String,
                )

                @test occursin(r"ONE", out) # Printed into STDOUT
                @test occursin(r"TWO", out) # Printed into STDOUT
                @test !occursin(r"RUN : 42", out) # Printed into the file
                @test filesize(one_out) != 0 # Created file
                @test filesize(two_out) != 0 # Created file
                @test occursin(r"RUN : 42", read(one_out, String)) # Printed into the file
                @test occursin(r"RUN : 42", read(two_out, String)) # Printed into the file

                _clean_up_outputs()
            end
        end
    end

    @testset "File -> File" begin

        @testset "--list & -a" begin
            _clean_up_outputs()
            out = read(`$julia $template $list_file --a 42 --list`, String)

            @test occursin(r"ONE", out) # Printed into STDOUT
            @test occursin(r"TWO", out) # Printed into STDOUT
            @test !occursin(r"RUN : 42", out) # Printed into the file
            @test filesize(one_out) != 0 # Created file
            @test filesize(two_out) != 0 # Created file
            @test occursin(r"RUN : 42", read(one_out, String)) # Printed into the file
            @test occursin(r"RUN : 42", read(two_out, String)) # Printed into the file

            _clean_up_outputs()
        end
    end

    @testset "-p (parallel)" begin

        @testset "File -> File (--list)" begin
            _clean_up_outputs()
            out = read(`$julia $template $list_file -p 2 --a 42 --list`, String)

            @test occursin(r"ONE", out) # Printed into STDOUT
            @test occursin(r"TWO", out) # Printed into STDOUT
            @test !occursin(r"RUN : 42", out) # Printed into the file
            @test filesize(one_out) != 0 # Created file
            @test filesize(two_out) != 0 # Created file
            @test occursin(r"RUN : 42", read(one_out, String)) # Printed into the file
            @test occursin(r"RUN : 42", read(two_out, String)) # Printed into the file

            _clean_up_outputs()
        end

        @static if Sys.isunix()

            @testset "STDIN -> File (--list)" begin
                _clean_up_outputs()
                out = read(
                    pipeline(`cat $list_file`, `$julia $template -p 2 --a 42 --list`),
                    String,
                )

                @test occursin(r"ONE", out) # Printed into STDOUT
                @test occursin(r"TWO", out) # Printed into STDOUT
                @test !occursin(r"RUN : 42", out) # Printed into the file
                @test filesize(one_out) != 0 # Created file
                @test filesize(two_out) != 0 # Created file
                @test occursin(r"RUN : 42", read(one_out, String)) # Printed into the file
                @test occursin(r"RUN : 42", read(two_out, String)) # Printed into the file
                _clean_up_outputs()
            end
        end
    end

    _clean_up()
end
