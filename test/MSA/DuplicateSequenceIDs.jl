using CodecZlib: GzipCompressor, transcode

@testset "Duplicate sequence identifiers" begin
    output_types = (
        AnnotatedMultipleSequenceAlignment,
        MultipleSequenceAlignment,
        NamedResidueMatrix{Matrix{Residue}},
        Matrix{Residue},
    )

    function records(format, ids)
        if format in (PIR, PIRSequences)
            join(">P1;$id\nTitle $i\nAC*\n" for (i, id) in enumerate(ids))
        elseif format === AnnotatedFASTASequences
            join(">$id\nAC\n01\n" for id in ids)
        else
            join(">$id\nAC\n" for id in ids)
        end
    end

    @testset "Record format $format" for format in (
        FASTA,
        A2M,
        A3M,
        PIR,
        FASTASequences,
        PIRSequences,
        AnnotatedFASTASequences,
    )
        outputs = format <: MSA.MSAFormat ? output_types : (nothing,)
        duplicated = records(format, ["dup", "other", "dup"])
        warning = "Sequence identifier dup is taken; using dup(1) instead."

        @testset "Output $T" for T in outputs
            args = T === nothing ? () : (T,)
            for input in (identity, IOBuffer)
                for options in ((;), (; fail_on_duplicate_seqnames = false))
                    parsed = @test_logs (:warn, warning) parse_file(
                        input(duplicated),
                        format,
                        args...;
                        options...,
                    )
                    if T === nothing
                        @test sequence_id.(parsed) == ["dup", "other", "dup(1)"]
                        @test getannotsequence(parsed[3], "OriginalSeqName") == "dup"
                    elseif T !== Matrix{Residue}
                        @test sequencenames(parsed) == ["dup", "other", "dup(1)"]
                        if T === AnnotatedMultipleSequenceAlignment
                            @test getannotsequence(parsed, "dup(1)", "OriginalSeqName") ==
                                  "dup"
                        end
                    else
                        @test size(parsed) == (3, 2)
                    end
                end

                duplicate_ids = [["dup", "other", "dup"], ["dup(1)", "dup(1)"]]
                # FastaIO rejects non-ASCII headers on streams before ID handling.
                if input === identity || !(format in (FASTA, A2M, A3M, FASTASequences))
                    push!(duplicate_ids, ["α", "β", "α"])
                end
                for ids in duplicate_ids
                    err = @test_logs try
                        parse_file(
                            input(records(format, ids)),
                            format,
                            args...;
                            fail_on_duplicate_seqnames = true,
                        )
                    catch e
                        e
                    end
                    @test err isa ArgumentError
                    @test sprint(showerror, err) ==
                          "ArgumentError: Duplicate sequence identifier: $(repr(last(ids)))."
                end

                # Literal suffixes are valid IDs, not evidence of a duplicate.
                unique = records(format, ["dup", "dup(1)", "dup(2)"])
                parsed = @test_logs parse_file(
                    input(unique),
                    format,
                    args...;
                    fail_on_duplicate_seqnames = true,
                )
                expected = parse_file(input(unique), format, args...)
                @test parsed == expected
                if T === AnnotatedMultipleSequenceAlignment
                    @test annotations(parsed) == annotations(expected)
                end
            end

            mktempdir() do dir
                for suffix in ("", ".gz")
                    path = joinpath(dir, "sequences" * suffix)
                    write(
                        path,
                        suffix == "" ? duplicated : transcode(GzipCompressor, duplicated),
                    )
                    @test_logs @test_throws ArgumentError read_file(
                        path,
                        format,
                        args...;
                        fail_on_duplicate_seqnames = true,
                    )
                    @test_logs (:warn, warning) read_file(path, format, args...)
                end
            end
        end

        # Exercise the default MSA output dispatch as well as explicit output types.
        @test_logs @test_throws ArgumentError parse_file(
            duplicated,
            format;
            fail_on_duplicate_seqnames = true,
        )
    end

    @testset "Sequence metadata and unequal lengths" begin
        fasta = ">one\nAC\n>two\nACDE\n"
        @test join.(parse_file(fasta, FASTASequences; fail_on_duplicate_seqnames = true)) ==
              ["AC", "ACDE"]
        for format in (PIRSequences, AnnotatedFASTASequences)
            parsed = parse_file(
                records(format, ["one", "two"]),
                format;
                fail_on_duplicate_seqnames = true,
            )
            if format === PIRSequences
                @test getannotsequence(parsed[2], "Title") == "Title 2"
                @test getannotsequence(parsed[2], "Type") == "P1"
            else
                @test getannotresidue(parsed[2], "Feature") == "01"
            end
        end
        msa = parse_file(
            ">one/2-3\nAC\n>two/5-6\nAC\n",
            FASTA;
            fail_on_duplicate_seqnames = true,
            generatemapping = true,
            useidcoordinates = true,
            keepinserts = true,
            deletefullgaps = false,
        )
        @test getsequencemapping(msa, "one/2-3") == [2, 3]
        @test getsequencemapping(msa, "two/5-6") == [5, 6]
    end

    @testset "Raw formats generate unique IDs" begin
        for args in ((Raw,), ((Raw, T) for T in output_types)..., (RawSequences,))
            @test parse_file("AC\nGT\n", args...; fail_on_duplicate_seqnames = true) ==
                  parse_file("AC\nGT\n", args...)
        end
    end

    @testset "Block format $format" for format in (Stockholm, Clustal)
        header = format === Stockholm ? "# STOCKHOLM 1.0\n" : "CLUSTAL\n\n"
        ending = format === Stockholm ? "//\n" : "\n"
        # Repeated annotation references are not sequence records.
        annotation = format === Stockholm ? "#=GR dup SS HH\n#=GS dup DE description\n" : ""
        block = "dup AC\n" * annotation * "other GT\n"
        valid = header * block * "\n" * block * ending
        duplicated = header * "dup AC\n" * annotation * "dup GT\n" * ending

        @testset "Output $T" for T in output_types
            for input in (identity, IOBuffer)
                @test_logs @test_throws ArgumentError parse_file(
                    input(duplicated),
                    format,
                    T;
                    fail_on_duplicate_seqnames = true,
                )
                # Default parsing still joins fragments carrying the same name.
                @test size(parse_file(input(duplicated), format, T)) == (1, 4)
                @test parse_file(
                    input(duplicated),
                    format,
                    T;
                    fail_on_duplicate_seqnames = false,
                ) == parse_file(input(duplicated), format, T)
                strict = @test_logs parse_file(
                    input(valid),
                    format,
                    T;
                    fail_on_duplicate_seqnames = true,
                )
                @test strict == parse_file(input(valid), format, T)
                @test size(strict) == (2, 4)
                if T === AnnotatedMultipleSequenceAlignment && format === Stockholm
                    @test getannotresidue(strict, "dup", "SS") == "HHHH"
                end
            end

            mktempdir() do dir
                for suffix in ("", ".gz")
                    path = joinpath(dir, "alignment" * suffix)
                    write_fixture(text) =
                        write(path, suffix == "" ? text : transcode(GzipCompressor, text))
                    write_fixture(duplicated)
                    @test_logs @test_throws ArgumentError read_file(
                        path,
                        format,
                        T;
                        fail_on_duplicate_seqnames = true,
                    )
                    write_fixture(valid)
                    @test read_file(path, format, T; fail_on_duplicate_seqnames = true) ==
                          parse_file(valid, format, T)

                    # Duplicate state must reset for each alignment in the same file.
                    write_fixture(valid * valid)
                    msas = eachmsa(path, format, T; fail_on_duplicate_seqnames = true)
                    @test length(@test_logs collect(msas)) == 2
                    @test !isopen(msas)

                    # Validation stays lazy and an error closes the iterator's stream.
                    write_fixture(valid * duplicated)
                    msas = eachmsa(path, format, T; fail_on_duplicate_seqnames = true)
                    @test first(msas) == parse_file(valid, format, T)
                    @test_logs @test_throws ArgumentError first(msas)
                    @test !isopen(msas)
                    @test !isopen(msas.io)
                end
            end
        end
    end

    @testset "Existing wrapped fixtures" begin
        for (format, file) in (
            (Stockholm, "clustalo-I20240512-trunc.aln-stockholm"),
            (Stockholm, "hmmer_multiblock.sto"),
            (Clustal, "PF09645.aln"),
            (Clustal, "PF09645.aln-num"),
        )
            path = joinpath(DATA, file)
            for T in output_types
                @test read_file(path, format, T; fail_on_duplicate_seqnames = true) ==
                      read_file(path, format, T)
            end
        end
    end

    @testset "Legacy custom loaders do not need the new keyword" begin
        struct LegacyDuplicateTestMSA <: MSA.MSAFormat end
        struct LegacyDuplicateTestSequences <: MSA.SequenceFormat end

        function MSA._load_sequences(
            io::Union{IO,AbstractString},
            ::Type{F};
            create_annotations::Bool = false,
        ) where {F<:Union{LegacyDuplicateTestMSA,LegacyDuplicateTestSequences}}
            ["one"], ["AC"], Annotations()
        end

        for options in ((;), (; fail_on_duplicate_seqnames = false))
            @test size(parse_file("", LegacyDuplicateTestMSA; options...)) == (1, 2)
            for T in output_types
                @test size(parse_file("", LegacyDuplicateTestMSA, T; options...)) == (1, 2)
            end
            @test join(only(parse_file("", LegacyDuplicateTestSequences; options...))) ==
                  "AC"
        end
    end
end
