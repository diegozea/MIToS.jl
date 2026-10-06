@testset "eachmsa" begin
    # Prepare plain and compressed fixtures without sharing helpers with other test files.
    _write_msa_fixture(path, contents) =
        write(path, endswith(path, ".gz") ? transcode(GzipCompressor, contents) : contents)

    @testset "Buffered line preservation" for compressed in (false, true)
        # Checking for another alignment must preserve long lines and special characters,
        # even with tiny buffers, repeated checks or an error.
        line = repeat("α", 10_000) * "\r\n"
        text = line * "tail"
        data = compressed ? transcode(GzipCompressor, text) : text
        wrapper = compressed ? Utils.GzipDecompressorStream : NoopStream
        io = wrapper(IOBuffer(data); bufsize = 1)
        @test_throws ArgumentError hasnextmsa(io, r"^other"; strict = true)
        @test hasnextmsa(io, r"^α"; strict = true)
        @test hasnextmsa(io, r"^α"; strict = true)
        first_line = MSA._read_msa_line(io)
        @test MSA._read_msa_line(io) == "tail"
        @test first_line == line
        @test eof(io)
        close(io)
    end

    @testset "$format" for format in (Stockholm, Clustal)
        @test support_eachmsa(format)
        multiple_warning = "Read only the first alignment; use `eachmsa` to read all."
        # Issue #202: different numbers of sequences and a repeated identifier across MSAs.
        first_record = """
        # STOCKHOLM 1.0
        #=GF SQ 4
        #=GF DE Example α
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
                # Read both alignments into each supported type. read_file must return
                # only the first alignment and warn that another follows.
                expected = [
                    parse_file(record, format, T; deletefullgaps = false) for
                    record in (first_record, second_record)
                ]
                for path in (plain, gzip)
                    msas = eachmsa(path, format, T; deletefullgaps = false)
                    @test eltype(msas) === T
                    alignments = @test_logs collect(msas)
                    @test alignments isa Vector{T}
                    @test size.(alignments) == [(4, 30), (3, 30)]
                    @test alignments == expected
                    msa = @test_logs (:warn, multiple_warning) read_file(path, format, T)
                    @test msa == parse_file(first_record, format, T)
                end
            end

            @testset "read_file warnings" begin
                expected = parse_file(first_record, format)
                header = format === Stockholm ? "# STOCKHOLM 1.0\n" : "CLUSTAL\n"
                for suffix in ("", ".gz")
                    path = joinpath(dir, "warning-test" * suffix)
                    # Trailing whitespace or comments must not cause a warning.
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

            @testset "Successive parsing, $(repr(newline))" for newline in ("\n", "\r\n")
                # Preserve columns to compare annotations without modification timestamps.
                options = (deletefullgaps = false,)
                expected = [
                    parse_file(record, format; options...) for
                    record in (first_record, second_record)
                ]
                # Reading one alignment from an open stream must leave the next readable.
                # Cover Unix and Windows line endings across several input streams,
                # with small buffers that split lines and headers across their boundaries.
                text = replace(contents, "\n" => newline)
                msa = @test_logs parse_file(text, format; options...)
                @test msa == first(expected)
                streams = (
                    IOBuffer(text),
                    NoopStream(IOBuffer(text); bufsize = 7),
                    Utils.GzipDecompressorStream(
                        IOBuffer(transcode(GzipCompressor, text));
                        bufsize = 7,
                    ),
                )
                if format === Clustal
                    _write_msa_fixture(gzip, text)
                    streams = (streams..., GZip.open(gzip))
                end
                for io in streams
                    try
                        msa = @test_logs parse_file(io, format; options...)
                        @test msa == first(expected)
                        @test annotations(msa) == annotations(first(expected))
                        if format === Stockholm || io isa TranscodingStream
                            @test position(io) ==
                                  sizeof(replace(first_record, "\n" => newline))
                        end
                        msa = @test_logs parse_file(io, format; options...)
                        @test msa == last(expected)
                    finally
                        close(io)
                    end
                end
            end

            @testset "Defaults, lookahead and cleanup" for path in (plain, gzip)
                # Defaults remove all-gap columns and return independently annotated MSAs.
                # Finishing must close the input while preserving the local source file.
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
                # Read alignments split across text blocks, passing options for mappings,
                # sequence coordinates, inserts and gap columns to each parser call.
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
                @testset "Pipe-backed IOStream" begin
                    # Clustal must read consecutive alignments even when the input cannot
                    # move backwards, with or without an additional buffer.
                    if Sys.isunix()
                        for buffered in (false, true)
                            fds = Vector{Cint}(undef, 2)
                            @test ccall(:pipe, Cint, (Ptr{Cint},), fds) == 0
                            io = Base.fdio(fds[1], true)
                            writer = Base.fdio(fds[2], true)
                            try
                                write(writer, contents)
                                close(writer)
                                @test_throws SystemError position(io)
                                input = buffered ? NoopStream(io; bufsize = 7) : io
                                @test parse_file(input, format) ==
                                      parse_file(first_record, format)
                                @test parse_file(input, format) ==
                                      parse_file(second_record, format)
                            finally
                                close(io)
                                close(writer)
                            end
                        end
                    end
                end
                @testset "Clustal headers, blocks and conservation" begin
                    # Handle header variants, residue counts and multiple text blocks.
                    # CLUSTAL_seq is a sequence name; conservation belongs to its own MSA.
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
                            second_msa = first(msas)
                            @test size(second_msa) == (1, 2)
                            @test stringsequence(second_msa, 1) == "KK"
                            @test getannotcolumn(second_msa, "cons", "") == ""
                            @test iterate(msas) === nothing
                        end
                    end
                end
                @testset "Conservation line endings" begin
                    # Preserve spaces in conservation annotations without including Unix
                    # or Windows newline characters when reading strings or files.
                    for newline in ("\n", "\r\n"), cons in ("*", "* ", " *")
                        record = join(
                            ("CLUSTAL", "", "a ACD", "b AEF", "  " * cons, ""),
                            newline,
                        )
                        @test getannotcolumn(parse_file(record, Clustal), "cons") == cons
                        for suffix in ("", ".gz")
                            path = joinpath(dir, "conservation.aln" * suffix)
                            _write_msa_fixture(path, record * record)
                            msa = @test_logs (:warn, multiple_warning) read_file(
                                path,
                                Clustal,
                            )
                            @test getannotcolumn(msa, "cons") == cons
                            @test getannotcolumn.(
                                collect(eachmsa(path, Clustal)),
                                "cons",
                            ) == [cons, cons]
                        end
                    end
                end
            end

            @testset "Empty, single and whitespace-separated records" begin
                # Accept empty files, extra whitespace and a missing final newline.
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
                # The do form must close the input after an early stop or an analysis error.
                # Explicitly closing a reader must also prevent further reading.
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
                # A malformed second alignment must not prevent reading the first.
                # Parsing errors close the input; invalid headers and missing files fail.
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
                # When the paired collection ends, the next alignment must remain unread.
                eachmsa(plain, format) do msas
                    @test length(collect(zip(msas, 1:1))) == 1
                    @test size(first(msas)) == (3, 30)
                end
            end

        end
    end

    @testset "User-defined format" begin
        # A format added through the public interface must work with gzip, parser options,
        # output types, warnings and cleanup after errors.
        # A test format giving the sequence count followed by that many sequence lines.
        struct _CountedMSAFormat <: MSAFormat end

        MSA.support_eachmsa(::Type{_CountedMSAFormat}) = true

        MSA.hasnextmsa(io::IO, ::Type{_CountedMSAFormat}; strict::Bool = false) =
            hasnextmsa(io, r"^\d+$"; strict = strict)

        function Utils.parse_file(
            io::IO,
            ::Type{_CountedMSAFormat},
            output::Type;
            kwargs...,
        )
            count = parse(Int, readline(io))
            sequences = join((readline(io) for _ = 1:count), '\n')
            parse_file(sequences, Raw, output; kwargs...)
        end

        @test support_eachmsa(_CountedMSAFormat)
        mktempdir() do dir
            record = "2\nA-C-\nATC-\n"
            options = (generatemapping = true, deletefullgaps = false)
            expected = parse_file("A-C-\nATC-", Raw; options...)
            for suffix in ("", ".gz")
                path = joinpath(dir, "counted.msas" * suffix)
                _write_msa_fixture(path, record * "\n1\nGG\n\n")
                eachmsa(path, _CountedMSAFormat; options...) do msas
                    @test !isempty(msas)
                    @test !isempty(msas)
                    msa = first(msas)
                    @test msa == expected
                    @test annotations(msa) == annotations(expected)
                    @test stringsequence(only(collect(msas)), 1) == "GG"
                    @test !isopen(msas)
                end
                # Output types and keywords go through the extension's public parser.
                matrices = collect(
                    eachmsa(
                        path,
                        _CountedMSAFormat,
                        Matrix{Residue};
                        deletefullgaps = false,
                    ),
                )
                @test size.(matrices) == [(2, 4), (1, 2)]
                msa = @test_logs (
                    :warn,
                    "Read only the first alignment; use `eachmsa` to read all.",
                ) read_file(path, _CountedMSAFormat; options...)
                @test annotations(msa) == annotations(expected)
                _write_msa_fixture(path, record * "bad header\n")
                msas = eachmsa(path, _CountedMSAFormat)
                @test size(first(msas)) == (2, 3)
                @test_throws ArgumentError iterate(msas)
                @test !isopen(msas)
                _write_msa_fixture(path, "\n\t\n")
                @test isempty(collect(eachmsa(path, _CountedMSAFormat)))
            end
        end
    end

    @testset "Non-iterable custom formats" begin
        # Formats without eachmsa support must retain their own read_file defaults.
        # A format with only a two-argument parser and no iterator support.
        struct _TwoArgumentMSAFormat <: MSAFormat end

        Utils.parse_file(io::IO, ::Type{_TwoArgumentMSAFormat}; kwargs...) =
            parse_file(io, Raw, MultipleSequenceAlignment; kwargs...)

        # A format whose parser defaults to an unannotated alignment.
        struct _DefaultOutputMSAFormat <: MSAFormat end

        # A lookahead method alone must not enable iteration or change read_file dispatch.
        MSA.hasnextmsa(io::IO, ::Type{_DefaultOutputMSAFormat}; strict::Bool = false) =
            error("This format has not opted into eachmsa.")

        Utils.parse_file(
            io::IO,
            ::Type{_DefaultOutputMSAFormat},
            output::Type = MultipleSequenceAlignment;
            kwargs...,
        ) = parse_file(io, Raw, output; kwargs...)

        mktempdir() do dir
            record = "A-C-\nATC-\n"
            options = (deletefullgaps = false,)
            expected = parse_file(record, Raw, MultipleSequenceAlignment; options...)
            for suffix in ("", ".gz")
                path = joinpath(dir, "custom.msa" * suffix)
                _write_msa_fixture(path, record)
                for format in (_TwoArgumentMSAFormat, _DefaultOutputMSAFormat)
                    @test !support_eachmsa(format)
                    @test_throws MethodError eachmsa(path, format)
                    msa = read_file(path, format; options...)
                    @test msa isa MultipleSequenceAlignment
                    @test msa == expected
                end
                msa = read_file(
                    path,
                    _DefaultOutputMSAFormat,
                    AnnotatedMultipleSequenceAlignment;
                    generatemapping = true,
                    options...,
                )
                @test msa isa AnnotatedMultipleSequenceAlignment
                @test getcolumnmapping(msa) == collect(1:4)
            end
        end
    end

    @testset "Unsupported format" begin
        # Both calling forms must reject an unsupported format before opening the file.
        @test !support_eachmsa(Raw)
        mktempdir() do dir
            path = joinpath(dir, "missing.txt")
            for args in ((path, Raw), (identity, path, Raw))
                @test_throws MethodError eachmsa(args...; deletefullgaps = false)
            end
        end
    end

    @testset "URL download lifetime" begin
        # Keep the downloaded file available while reading, then delete it when the reader
        # reaches the end or is explicitly closed.
        base = "https://raw.githubusercontent.com/diegozea/MIToS.jl/0ce717038b642d550f710ba7ea095d791812ff6e/"
        @test_throws MethodError eachmsa(base * "test/data/PF09645_full.stockholm", Raw)
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
