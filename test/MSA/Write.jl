@testset "Writing multiple alignments" begin
    # An alignment wrapper with custom I/O methods that also iterates over its residues.
    struct _CustomMSA <: AbstractMatrix{Residue}
        msa::AnnotatedMultipleSequenceAlignment
    end

    Base.size(custom::_CustomMSA) = size(custom.msa)
    Base.getindex(custom::_CustomMSA, i::Int, j::Int) = custom.msa[i, j]

    function Utils.parse_file(
        io::IO,
        format::Type{<:Union{Stockholm,Clustal}},
        ::Type{_CustomMSA};
        kwargs...,
    )
        header = format === Stockholm ? "# STOCKHOLM" : "CLUSTAL"
        startswith(readline(io), header) || error("Missing alignment header.")
        _CustomMSA(parse_file(io, format; kwargs...))
    end

    Utils.print_file(io::IO, custom::_CustomMSA, ::Type{Stockholm}) =
        print_file(io, custom.msa, Stockholm)

    # Also cover a writer specialized on the actual file stream type.
    Utils.print_file(
        io::Union{IOStream,Utils.GzipCompressorStream},
        custom::_CustomMSA,
        ::Type{Clustal},
    ) = print_file(io, custom.msa, Clustal)

    # A non-matrix object with a writer, excluded from alignment collections.
    struct _NonMatrixMSA end

    Utils.print_file(io::IO, ::_NonMatrixMSA, ::Type{<:Union{Stockholm,Clustal}}) = nothing

    records = (
        "# STOCKHOLM 1.0\n#=GF ID first\n#=GS a DE first sequence\na AC-\nb ADE\n#=GC cons * .\n//\n",
        "# STOCKHOLM 1.0\n#=GF ID second\na WK\n#=GC cons :*\n//\n",
    )
    msas = [parse_file(record, Stockholm; deletefullgaps = false) for record in records]
    msa_types = (
        Matrix{Residue},
        NamedResidueMatrix{Matrix{Residue}},
        MultipleSequenceAlignment,
        AnnotatedMultipleSequenceAlignment,
    )

    @testset "$format" for format in (Stockholm, Clustal)
        mktempdir() do dir
            @testset "Custom I/O methods" for suffix in ("", ".gz")
                source = joinpath(dir, "source" * suffix)
                output = joinpath(dir, "custom" * suffix)
                write_file(source, msas[1], format)
                expected = read_file(source, format; generatemapping = true)
                write_file(output, _CustomMSA(msas[1]), format)
                custom = read_file(output, format, _CustomMSA; generatemapping = true)
                @test custom.msa == expected
                @test annotations(custom.msa) == annotations(expected)
                write_file(source, _CustomMSA.(msas), format)
                @test_logs (
                    :warn,
                    "Read only the first alignment; use `eachmsa` to read all.",
                ) read_file(source, format, _CustomMSA; generatemapping = true)
                eachmsa(source, format, _CustomMSA; generatemapping = true) do input
                    write_file(output, input, format)
                end
                actual = collect(eachmsa(output, format, _CustomMSA))
                @test getfield.(actual, :msa) == msas
                if format === Stockholm
                    @test annotations(actual[1].msa) == annotations(expected)
                end
                mixed = Union{_CustomMSA,AnnotatedMultipleSequenceAlignment}[
                    _CustomMSA(msas[1]),
                    msas[2],
                ]
                write_file(output, mixed, format)
                @test collect(eachmsa(output, format)) == msas
            end

            @testset "Alignment types" for T in msa_types
                typed_msas = [
                    parse_file(record, Stockholm, T; deletefullgaps = false) for
                    record in records
                ]
                expected = join(sprint(print_file, msa, format) for msa in typed_msas)
                @test sprint(print_file, typed_msas, format) == expected
                for suffix in ("", ".gz")
                    path = joinpath(dir, "alignments" * suffix)
                    @test write_file(path, typed_msas, format) === nothing
                    actual = collect(eachmsa(path, format, T; deletefullgaps = false))
                    @test actual == typed_msas
                    @test sequencenames.(actual) == sequencenames.(typed_msas)
                    if T === AnnotatedMultipleSequenceAlignment
                        if format === Stockholm
                            @test annotations.(actual) == annotations.(typed_msas)
                        else
                            @test getannotcolumn.(actual, "cons") ==
                                  getannotcolumn.(typed_msas, "cons")
                        end
                    end
                end
            end

            @testset "eachmsa input" for suffix in ("", ".gz")
                source = joinpath(dir, "input" * suffix)
                output = joinpath(dir, "output" * suffix)
                write_file(source, msas, format)
                eachmsa(source, format) do input
                    write_file(output, input, format)
                end
                @test collect(eachmsa(output, format)) == msas
            end

            @testset "Empty input" for empty in
                                       (AnnotatedMultipleSequenceAlignment[], (), Union{}[])
                @test sprint(print_file, empty, format) == ""
                for suffix in ("", ".gz")
                    output = joinpath(dir, "empty" * suffix)
                    write_file(output, empty, format)
                    @test isempty(collect(eachmsa(output, format)))
                end
            end

            @testset "Append" begin
                path = joinpath(dir, "append")
                write_file(path, msas[1], format)
                write_file(path, msas[2:2], format, "a")
                @test collect(eachmsa(path, format)) == msas
            end

            @testset "Incremental output" begin
                io = IOBuffer()
                steps = Int[]
                input = Iterators.filter(msas) do msa
                    push!(steps, length(steps) + 1)
                    if length(steps) == 2
                        @test String(take!(io)) == sprint(print_file, msas[1], format)
                    end
                    true
                end
                @test print_file(io, input, format) === nothing
                @test steps == [1, 2]
                @test String(take!(io)) == sprint(print_file, msas[2], format)
                @test isopen(io)
            end

            @testset "Iteration errors close output" for suffix in ("", ".gz")
                path = joinpath(dir, "partial" * suffix)
                input = Iterators.filter(msas) do msa
                    msa === first(msas) || error("iteration failed")
                end
                @test_throws ErrorException("iteration failed") write_file(
                    path,
                    input,
                    format,
                )
                @test collect(eachmsa(path, format)) == msas[1:1]
            end

            @testset "Reject missing element types" begin
                invalid_inputs = (
                    (Any[msas[1]], "Expected an alignment `eltype`; got `Any`."),
                    (Any[], "Expected an alignment `eltype`; got `Any`."),
                    (
                        (error("must not iterate") for _ = 1:1),
                        "Iterator must declare an alignment `eltype`.",
                    ),
                )
                for (input, message) in invalid_inputs
                    path = joinpath(dir, "invalid")
                    @test_throws ArgumentError(message) write_file(path, input, format)
                    io = IOBuffer()
                    @test_throws ArgumentError(message) print_file(io, input, format)
                    @test position(io) == 0
                end
            end

            @testset "Reject non-alignment elements" for suffix in ("", ".gz")
                path = joinpath(dir, "unsupported" * suffix)
                recursive = Vector[]
                push!(recursive, recursive)
                for input in (
                    1,
                    'A',
                    [1],
                    ['A'],
                    recursive,
                    [msas],
                    [ones(Int, 1, 1)],
                    [_NonMatrixMSA()],
                )
                    @test_throws MethodError write_file(path, input, format)
                end
                mixed = Union{AnnotatedMultipleSequenceAlignment,Int}[msas[1], 1]
                @test_throws MethodError write_file(path, mixed, format)
                @test read_file(path, format) == msas[1]
            end
        end
    end

    @testset "Pfam example preserves mappings" begin
        fixture = joinpath(DATA, "PF09645_full.stockholm")
        expected = read_file(fixture, Stockholm)
        mktempdir() do dir
            source = joinpath(dir, "msas.sto.gz")
            output = joinpath(dir, "aligned-msas.sto.gz")
            write(
                source,
                transcode(GzipCompressor, repeat(read(fixture, String) * "\n", 2)),
            )
            eachmsa(
                source,
                Stockholm;
                generatemapping = true,
                useidcoordinates = true,
            ) do input
                write_file(output, input, Stockholm)
            end
            saved = @test_logs collect(eachmsa(output, Stockholm))
            @test saved == [expected, expected]
            for msa in saved
                @test getcolumnmapping(msa) == 6:115
                @test getannotfile(msa, "NCol") == "120"
                @test getsequencemapping(msa, 1)[1:4] == [0, 0, 0, 3]
                @test getsequencemapping(msa, "F112_SSV1/3-112") == 3:112
                @test getannotcolumn(msa) == getannotcolumn(expected)
                @test getannotresidue(msa) == getannotresidue(expected)
                @test getannotsequence(msa, "F112_SSV1/3-112", "DR") == "PDB; 2VQC A; 4-73;"
            end
            eachmsa(output, Stockholm) do input
                write_file(source, input, Stockholm)
            end
            reread = @test_logs collect(eachmsa(source, Stockholm))
            @test reread == saved
            @test annotations.(reread) == annotations.(saved)
        end
    end

    @testset "Format dispatch and keywords" begin
        @test !applicable(print_file, IOBuffer(), msas, FASTA)
        @test sprint(io -> print_file(io, msas, Clustal; showcounts = true)) == join(
            sprint(io -> print_file(io, msa, Clustal; showcounts = true)) for msa in msas
        )
    end
end
