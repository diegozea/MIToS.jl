using CodecZlib: GzipCompressor, transcode
import Downloads

"""
Write an alignment fixture, compressing it when the filename ends in `.gz`.
"""
_write_msa_fixture(path, contents) =
    write(path, endswith(path, ".gz") ? transcode(GzipCompressor, contents) : contents)

@testset "eachmsa" begin
    @testset "$format" for format in (Stockholm, Clustal)
        multiple_warning = "Read only the first alignment; use `eachmsa` to read all."
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
            _write_msa_fixture(plain, contents)
            _write_msa_fixture(gzip, contents)

            @testset "Plain and gzip, output $T" for T in output_types
                for path in (plain, gzip)
                    msas = eachmsa(path, format, T; deletefullgaps = false)
                    @test eltype(msas) === T
                    alignments = @test_logs collect(msas)
                    @test alignments isa Vector{T}
                    @test size.(alignments) == [(4, 30), (3, 30)]
                    @test alignments[1] ==
                          parse_file(first_record, format, T; deletefullgaps = false)
                    @test alignments[2] ==
                          parse_file(second_record, format, T; deletefullgaps = false)
                end
            end

            @testset "read_file warnings" begin
                for T in output_types, path in (plain, gzip)
                    msa = @test_logs (:warn, multiple_warning) read_file(path, format, T)
                    @test msa == parse_file(first_record, format, T)
                end
                expected = parse_file(first_record, format)
                header = format === Stockholm ? "# STOCKHOLM 1.0\n" : "CLUSTAL\n"
                for suffix in ("", ".gz")
                    path = joinpath(dir, "warning-test" * suffix)
                    for trailing in ("", "\n \t\n", "\n# end of file\n")
                        text = first_record * trailing
                        _write_msa_fixture(path, text)
                        msa = @test_logs read_file(path, format)
                        @test msa == expected
                    end
                    # Warn once, even with three alignments. Detecting another header
                    # must not require reading or validating its alignment.
                    for following in (
                        second_record * second_record,
                        header,
                        header * "not a valid alignment\n",
                    )
                        text = first_record * "\n \t\n" * following
                        _write_msa_fixture(path, text)
                        msa = @test_logs (:warn, multiple_warning) read_file(path, format)
                        @test msa == expected
                    end
                end
            end

            @testset "Successive parse_file calls" begin
                expected = parse_file(first_record, format)
                msa = @test_logs parse_file(contents, format)
                @test msa == expected
                # Reading one alignment from an open stream must leave the next readable.
                io = IOBuffer(contents)
                msa = @test_logs parse_file(io, format)
                @test msa == expected
                if format === Stockholm
                    @test position(io) == sizeof(first_record)
                end
                msa = @test_logs parse_file(io, format)
                @test msa == parse_file(second_record, format)
            end

            @testset "Defaults, lookahead and cleanup" for path in (plain, gzip)
                msas = eachmsa(path, format)
                @test isopen(msas)
                @test position(msas.io) == 0 # opening does not parse the first MSA
                @test !isempty(msas)
                @test !isempty(msas) # repeated lookahead must not lose a record
                alignments = collect(msas)
                @test alignments isa Vector{AnnotatedMultipleSequenceAlignment}
                @test size.(alignments) == [(4, 29), (3, 30)]
                if format === Stockholm
                    @test getannotfile.(alignments, "SQ") == ["4", "3"]
                end
                @test stringsequence(alignments[2], 1) == "LPENWQALLDDTGTYFYANHLTKTSQWEHP"
                setannotfile!(alignments[1], "ID", "changed")
                @test getannotfile(alignments[2], "ID", "") == ""
                @test !isopen(msas)
                @test !isopen(msas.io)
                @test isempty(collect(msas))
                @test close(msas) === nothing
                @test isfile(path) # never remove a user's local source
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
                        _write_msa_fixture(path, payload)
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
                    _write_msa_fixture(path * suffix, payload)
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
            end

            @testset "zip does not consume an extra MSA" begin
                eachmsa(plain, format) do msas
                    @test length(collect(zip(msas, 1:1))) == 1
                    @test size(first(msas)) == (3, 30)
                end
            end

        end
    end

    @testset "Unsupported format" begin
        mktempdir() do dir
            @test_throws ArgumentError eachmsa(joinpath(dir, "missing.txt"), Raw)
        end
    end

    @testset "URL download lifetime" begin
        base = "https://raw.githubusercontent.com/diegozea/MIToS.jl/0ce717038b642d550f710ba7ea095d791812ff6e/"
        @test_throws ArgumentError eachmsa(base * "test/data/PF09645_full.stockholm", Raw)
        for (format, file) in (
            (Stockholm, "test/data/PF09645_full.stockholm"),
            (Stockholm, "docs/data/PF18883.stockholm.gz"),
            (Clustal, "test/data/PF09645.aln"),
        )
            url = base * file
            for exhaust in (false, true)
                msas = try
                    eachmsa(url, format)
                catch err
                    # Skip only proxy/DNS, connection, or timeout failures.
                    err isa Downloads.RequestError && err.code in (5, 6, 7, 28) ||
                        rethrow()
                    @test_skip eachmsa(url, format)
                    break
                end
                temporary = msas.temporary
                try
                    @test isfile(temporary)
                    @test first(msas) == read_file(joinpath(DATA, "..", "..", file), format)
                    @test isfile(temporary)
                    if exhaust
                        @test isempty(msas)
                    else
                        close(msas)
                    end
                    @test !isfile(temporary)
                    @test !isopen(msas.io)
                finally
                    close(msas)
                end
            end
        end
    end
end
