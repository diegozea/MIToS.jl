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
    @testset "$format" for format in (Stockholm, Clustal)
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
        if format === Clustal
            first_record, second_record = map((first_record, second_record)) do record
                io = IOBuffer()
                print_file(
                    io,
                    parse_file(
                        record,
                        Stockholm,
                        MultipleSequenceAlignment;
                        deletefullgaps = false,
                    ),
                    Clustal,
                )
                String(take!(io))
            end
        end
        contents = first_record * second_record
        output_types = (
            AnnotatedMultipleSequenceAlignment,
            MultipleSequenceAlignment,
            NamedResidueMatrix{Matrix{Residue}},
            Matrix{Residue},
        )
        mktempdir() do dir
            plain = joinpath(dir, format === Stockholm ? "multiple.sto" : "multiple.aln")
            gzip = plain * ".gz"
            write(plain, contents)
            write(gzip, transcode(GzipCompressor, contents))

            @testset "Plain and gzip, output $T" for T in output_types
                for path in (plain, gzip)
                    msas = eachmsa(path, format, T; deletefullgaps = false)
                    @test eltype(msas) === T
                    @test isopen(msas)
                    @test position(msas.io) == 0 # opening does not parse the first MSA
                    @test !isempty(msas)
                    @test !isempty(msas) # repeated lookahead must not lose a record
                    alignments = collect(msas)
                    @test alignments isa Vector{T}
                    @test size.(alignments) == [(4, 30), (3, 30)]
                    @test alignments[1] ==
                          parse_file(first_record, format, T; deletefullgaps = false)
                    @test alignments[2] ==
                          parse_file(second_record, format, T; deletefullgaps = false)
                    @test !isopen(msas)
                    @test !isopen(msas.io)
                    @test isempty(collect(msas))
                    @test close(msas) === nothing
                    @test isfile(path) # never remove a user's local source
                    if format === Stockholm
                        @test read_file(path, format, T) ==
                              parse_file(first_record, format, T)
                    end
                end
            end

            @testset "Defaults and independent annotations" begin
                alignments = collect(eachmsa(plain, format))
                @test alignments isa Vector{AnnotatedMultipleSequenceAlignment}
                @test size.(alignments) == [(4, 29), (3, 30)]
                if format === Stockholm
                    @test getannotfile.(alignments, "SQ") == ["4", "3"]
                end
                @test stringsequence(alignments[2], 1) == "LPENWQALLDDTGTYFYANHLTKTSQWEHP"
                setannotfile!(alignments[1], "ID", "changed")
                @test getannotfile(alignments[2], "ID", "") == ""
            end

            @testset "Wrapped records and parsing keywords" begin
                names =
                    format === Stockholm ?
                    ("PF09645_full.stockholm", "hmmer_multiblock.sto") :
                    ("PF09645.aln", "PF09645.aln-num")
                record, wrapped = map(name -> read(joinpath(DATA, name), String), names)
                path = joinpath(dir, "wrapped.sto")
                write(path, wrapped * "\n" * record)
                options = (
                    generatemapping = true,
                    useidcoordinates = true,
                    keepinserts = true,
                    deletefullgaps = false,
                )
                alignments = collect(eachmsa(path, format; options...))
                @test length(alignments) == 2
                for (msa, text) in zip(alignments, (wrapped, record))
                    expected = parse_file(text, format; options...)
                    @test msa == expected
                    @test annotations(msa) == annotations(expected)
                    @test getcolumnmapping(msa) == getcolumnmapping(expected)
                end
            end

            if format === Clustal
                @testset "Clustal headers, blocks and conservation" begin
                    wrapped = """
                    CLUSTAL W (2.1) multiple sequence alignment

                    CLUSTAL_seq A-C 2
                    b           A-D 2
                                * .

                    CLUSTAL_seq EFG 5
                    b           EYG 5
                                * *

                    """
                    next_record = "CLUSTALW (1.83) multiple sequence alignment\n\nCLUSTAL_seq KK 2\n"
                    for suffix in ("", ".gz")
                        path = joinpath(dir, "headers.aln" * suffix)
                        payload = wrapped * next_record
                        write(
                            path,
                            isempty(suffix) ? payload : transcode(GzipCompressor, payload),
                        )
                        eachmsa(path, Clustal; deletefullgaps = false) do msas
                            first_msa = first(msas)
                            @test size(first_msa) == (2, 6)
                            @test stringsequence(first_msa, 1) == "A-CEFG"
                            @test getannotcolumn(first_msa, "cons") == "* .* *"
                            @test !isempty(msas)
                            @test !isempty(msas)
                            second_msa = first(msas)
                            @test size(second_msa) == (1, 2)
                            @test stringsequence(second_msa, 1) == "KK"
                            @test getannotcolumn(second_msa, "cons", "") == ""
                            @test iterate(msas) === nothing
                        end
                    end
                end
            end

            @testset "Empty, single and whitespace-separated records" begin
                path = joinpath(dir, "edge.sto")
                for text in ("", " \n\t\r\n")
                    write(path, text)
                    @test isempty(collect(eachmsa(path, format)))
                end
                write(path, chomp(first_record)) # no final newline
                @test length(collect(eachmsa(path, format))) == 1
                write(path, "\n \t\n" * first_record * "\r\n\t\n" * second_record * " \t\n")
                @test size.(collect(eachmsa(path, format))) == [(4, 29), (3, 30)]
            end

            @testset "Early termination and user exceptions" begin
                for path in (plain, gzip)
                    reader = nothing
                    result = eachmsa(path, format) do msas
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
                    @test_throws ErrorException eachmsa(path, format) do msas
                        reader = msas
                        first(msas)
                        error("analysis failed")
                    end
                    @test !isopen(reader.io)
                    reader = eachmsa(path, format)
                    @test size(first(reader)) == (4, 29)
                    close(reader)
                    @test !isopen(reader.io)
                    @test iterate(reader) === nothing
                end
            end

            @testset "Lazy errors and cleanup" begin
                path = joinpath(dir, "malformed.sto")
                malformed =
                    format === Stockholm ? "# STOCKHOLM 1.0\na AAA\nb A\n//\n" :
                    "CLUSTAL\n\na AAA\nb A\n"
                for suffix in ("", ".gz")
                    payload = first_record * malformed
                    write(
                        path * suffix,
                        isempty(suffix) ? payload : transcode(GzipCompressor, payload),
                    )
                    msas = eachmsa(path * suffix, format)
                    @test size(first(msas)) == (4, 29)
                    @test_throws ErrorException iterate(msas)
                    @test !isopen(msas.io)
                end
                write(path, "invalid header\n" * second_record)
                msas = eachmsa(path, format)
                @test_throws ArgumentError iterate(msas)
                @test !isopen(msas.io)
                if format === Stockholm
                    write(path, first_record * "invalid header\n")
                    msas = eachmsa(path, format)
                    @test size(first(msas)) == (4, 29)
                    @test_throws ArgumentError iterate(msas)
                    @test !isopen(msas.io)
                end
                @test_throws SystemError eachmsa(joinpath(dir, "missing.sto"), format)
                @test_throws ArgumentError eachmsa(plain, Raw)
                @test_throws ArgumentError eachmsa(joinpath(dir, "missing.txt"), Raw)
            end

            @testset "zip does not consume an extra MSA" begin
                eachmsa(plain, format) do msas
                    @test length(collect(zip(msas, 1:1))) == 1
                    @test size(first(msas)) == (3, 30)
                end
            end

            @testset "URL download lifetime" begin
                for suffix in ("", ".gz")
                    payload = read(isempty(suffix) ? plain : gzip)
                    _with_msa_server(payload) do base, requests
                        url = base * "/multiple.sto" * suffix
                        @test_throws ArgumentError eachmsa(url, Raw)
                        @test isempty(requests) # unsupported formats must not download the URL
                        msas = eachmsa(url, format)
                        temporary = msas.temporary
                        @test isfile(temporary)
                        @test length(requests) == 1
                        @test size(first(msas)) == (4, 29)
                        @test isfile(temporary)
                        @test size.(collect(msas)) == [(3, 30)]
                        @test length(requests) == 1 # no redownload between records
                        @test !isfile(temporary)
                        @test !isopen(msas.io)
                        eachmsa(url, format) do reader
                            temporary = reader.temporary
                            first(reader)
                        end
                        @test !isfile(temporary)
                        @test length(requests) == 2
                    end
                end
                _with_msa_server(
                    Vector{UInt8}(codeunits("invalid header\n")),
                ) do base, requests
                    msas = eachmsa(base * "/bad.sto", format)
                    temporary = msas.temporary
                    @test_throws ArgumentError first(msas)
                    @test !isfile(temporary)
                    @test !isopen(msas.io)
                end
            end
        end
    end
end
