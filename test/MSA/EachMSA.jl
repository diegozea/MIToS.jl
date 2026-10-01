using CodecZlib: GzipCompressor, transcode
using Sockets

# Serve an in-memory fixture locally so URL tests also count downloads and do not
# depend on an external service. Each request uses a fresh connection.
function _with_msa_server(f, payload)
    server = listen(ip"127.0.0.1", 0)
    port = getsockname(server)[2]
    requests = String[]
    task = @async begin
        while isopen(server)
            socket = try
                accept(server)
            catch
                isopen(server) && rethrow()
                break
            end
            try
                push!(requests, readline(socket))
                while !isempty(readline(socket))
                end
                write(
                    socket,
                    "HTTP/1.1 200 OK\r\nContent-Length: $(length(payload))\r\nConnection: close\r\n\r\n",
                )
                write(socket, payload)
                flush(socket)
            finally
                close(socket)
            end
        end
    end
    try
        return f("http://127.0.0.1:$port", requests)
    finally
        close(server)
        wait(task)
    end
end

@testset "eachmsa" begin
    # Issue #202: different numbers of sequences and a repeated identifier across MSAs.
    first_record = """
    # STOCKHOLM 1.0
    #=GF SQ 4
    0.0 --EKWVKSKDEEGAYYYHDQGTNEVRWEKP
    1.7656331811711743 ---GWTVFRTKNGAAYYVETRTQQATWENP
    3.5312663623423486 -HPPWVQRMTPAGRTYYHYLQTRETTWTDP
    5.296899543513523 -AS-WEAAKDHDGHTYYIHEIKRETSWEIP
    //
    """
    second_record = """
    # STOCKHOLM 1.0
    #=GF SQ 3
    0.0 LPENWQALLDDTGTYFYANHLTKTSQWEHP
    0.1088478902405901 LPEVWQALLDDTGTVFYINHLTKTSQWEHP
    0.2176957804811802 LPEEWQALLDGTGTVFYINHGTKTSQWEHP
    //
    """
    contents = first_record * second_record
    output_types = (
        AnnotatedMultipleSequenceAlignment,
        MultipleSequenceAlignment,
        NamedResidueMatrix{Matrix{Residue}},
        Matrix{Residue},
    )
    mktempdir() do dir
        plain = joinpath(dir, "multiple.sto")
        gzip = plain * ".gz"
        write(plain, contents)
        write(gzip, transcode(GzipCompressor, contents))

        @testset "Plain and gzip, output $T" for T in output_types
            for path in (plain, gzip)
                msas = eachmsa(path, Stockholm, T; deletefullgaps = false)
                @test eltype(msas) === T
                @test isopen(msas)
                @test position(msas.io) == 0 # opening does not parse the first MSA
                @test !isempty(msas)
                @test !isempty(msas) # repeated lookahead must not lose a record
                alignments = collect(msas)
                @test alignments isa Vector{T}
                @test size.(alignments) == [(4, 30), (3, 30)]
                @test alignments[1] ==
                      parse_file(first_record, Stockholm, T; deletefullgaps = false)
                @test alignments[2] ==
                      parse_file(second_record, Stockholm, T; deletefullgaps = false)
                @test !isopen(msas)
                @test !isopen(msas.io)
                @test isempty(collect(msas))
                @test close(msas) === nothing
                @test isfile(path) # never remove a user's local source
                @test read_file(path, Stockholm, T) ==
                      parse_file(first_record, Stockholm, T)
            end
        end

        @testset "Defaults and independent annotations" begin
            alignments = collect(eachmsa(plain, Stockholm))
            @test alignments isa Vector{AnnotatedMultipleSequenceAlignment}
            @test size.(alignments) == [(4, 29), (3, 30)]
            @test getannotfile.(alignments, "SQ") == ["4", "3"]
            @test stringsequence(alignments[2], 1) == "LPENWQALLDDTGTYFYANHLTKTSQWEHP"
            setannotfile!(alignments[1], "ID", "changed")
            @test getannotfile(alignments[2], "ID", "") == ""
        end

        @testset "Wrapped records and parsing keywords" begin
            record = read(joinpath(DATA, "PF09645_full.stockholm"), String)
            wrapped = read(joinpath(DATA, "hmmer_multiblock.sto"), String)
            path = joinpath(dir, "wrapped.sto")
            write(path, wrapped * "\n" * record)
            options = (
                generatemapping = true,
                useidcoordinates = true,
                keepinserts = true,
                deletefullgaps = false,
            )
            alignments = collect(eachmsa(path, Stockholm; options...))
            @test length(alignments) == 2
            for (msa, text) in zip(alignments, (wrapped, record))
                expected = parse_file(text, Stockholm; options...)
                @test msa == expected
                @test annotations(msa) == annotations(expected)
                @test getcolumnmapping(msa) == getcolumnmapping(expected)
            end
        end

        @testset "Empty, single and whitespace-separated records" begin
            path = joinpath(dir, "edge.sto")
            for text in ("", " \n\t\r\n")
                write(path, text)
                @test isempty(collect(eachmsa(path, Stockholm)))
            end
            write(path, chomp(first_record)) # no final newline
            @test length(collect(eachmsa(path, Stockholm))) == 1
            write(path, "\n \t\n" * first_record * "\r\n\t\n" * second_record * " \t\n")
            @test size.(collect(eachmsa(path, Stockholm))) == [(4, 29), (3, 30)]
        end

        @testset "Early termination and user exceptions" begin
            for path in (plain, gzip)
                reader = nothing
                result = eachmsa(path, Stockholm) do msas
                    reader = msas
                    for msa in msas
                        @test size(msa) == (4, 29)
                        break
                    end
                    :finished
                end
                @test result === :finished
                @test !isopen(reader.io)
                @test iterate(reader) === nothing
                @test_throws ErrorException eachmsa(path, Stockholm) do msas
                    reader = msas
                    first(msas)
                    error("analysis failed")
                end
                @test !isopen(reader.io)
                reader = eachmsa(path, Stockholm)
                @test size(first(reader)) == (4, 29)
                close(reader)
                @test !isopen(reader.io)
                @test iterate(reader) === nothing
            end
        end

        @testset "Lazy errors and cleanup" begin
            path = joinpath(dir, "malformed.sto")
            malformed = "# STOCKHOLM 1.0\na AAA\nb A\n//\n"
            for suffix in ("", ".gz")
                payload = first_record * malformed
                write(
                    path * suffix,
                    isempty(suffix) ? payload : transcode(GzipCompressor, payload),
                )
                msas = eachmsa(path * suffix, Stockholm)
                @test size(first(msas)) == (4, 29)
                @test_throws ErrorException iterate(msas)
                @test !isopen(msas.io)
            end
            write(path, first_record * "not a Stockholm header\n")
            msas = eachmsa(path, Stockholm)
            @test size(first(msas)) == (4, 29)
            @test_throws ArgumentError iterate(msas)
            @test !isopen(msas.io)
            @test_throws SystemError eachmsa(joinpath(dir, "missing.sto"), Stockholm)
            @test_throws MethodError eachmsa(plain, Clustal)
        end

        @testset "zip does not consume an extra MSA" begin
            eachmsa(plain, Stockholm) do msas
                @test length(collect(zip(msas, 1:1))) == 1
                @test size(first(msas)) == (3, 30)
            end
        end

        @testset "URL download lifetime" begin
            for suffix in ("", ".gz")
                payload = read(isempty(suffix) ? plain : gzip)
                _with_msa_server(payload) do base, requests
                    url = base * "/multiple.sto" * suffix
                    msas = eachmsa(url, Stockholm)
                    temporary = msas.temporary
                    @test isfile(temporary)
                    @test length(requests) == 1
                    @test size(first(msas)) == (4, 29)
                    @test isfile(temporary)
                    @test size.(collect(msas)) == [(3, 30)]
                    @test length(requests) == 1 # no redownload between records
                    @test !isfile(temporary)
                    @test !isopen(msas.io)
                    eachmsa(url, Stockholm) do reader
                        temporary = reader.temporary
                        first(reader)
                    end
                    @test !isfile(temporary)
                    @test length(requests) == 2
                end
            end
            _with_msa_server(Vector{UInt8}(codeunits("invalid header\n"))) do base, requests
                msas = eachmsa(base * "/bad.sto", Stockholm)
                temporary = msas.temporary
                @test_throws ArgumentError first(msas)
                @test !isfile(temporary)
                @test !isopen(msas.io)
            end
        end
    end
end
