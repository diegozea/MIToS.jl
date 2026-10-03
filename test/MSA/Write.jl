"""
A typed iterator without a length, calling `before_next` before each iteration step.
"""
struct _MSAWriteIterator{T,F}
    msas::Vector{T}
    before_next::F
end

Base.IteratorSize(::Type{<:_MSAWriteIterator}) = Base.SizeUnknown()
Base.eltype(::Type{<:_MSAWriteIterator{T}}) where {T} = T

function Base.iterate(msas::_MSAWriteIterator, state::Int = 1)
    msas.before_next(state)
    iterate(msas.msas, state)
end

@testset "Writing multiple alignments" begin
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

            @testset "eachmsa and empty input" for suffix in ("", ".gz")
                source = joinpath(dir, "input" * suffix)
                output = joinpath(dir, "output" * suffix)
                write_file(source, msas, format)
                eachmsa(source, format) do input
                    write_file(output, input, format)
                end
                @test collect(eachmsa(output, format)) == msas
                write_file(output, AnnotatedMultipleSequenceAlignment[], format)
                @test isempty(collect(eachmsa(output, format)))
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
                input = _MSAWriteIterator(
                    msas,
                    state -> begin
                        push!(steps, state)
                        if state == 2
                            @test String(take!(io)) == sprint(print_file, msas[1], format)
                        end
                    end,
                )
                @test print_file(io, input, format) === nothing
                @test steps == [1, 2, 3]
                @test String(take!(io)) == sprint(print_file, msas[2], format)
                @test isopen(io)
            end

            @testset "Iteration errors close output" for suffix in ("", ".gz")
                path = joinpath(dir, "partial" * suffix)
                input = _MSAWriteIterator(
                    msas,
                    state -> begin
                        state == 2 && error("iteration failed")
                    end,
                )
                @test_throws ErrorException("iteration failed") write_file(
                    path,
                    input,
                    format,
                )
                @test collect(eachmsa(path, format)) == msas[1:1]
            end

            @testset "Reject element types before writing" begin
                invalid_inputs = (
                    (Any[msas[1]], "Expected an alignment eltype; got Any."),
                    (Any[], "Expected an alignment eltype; got Any."),
                    ([1], "Expected an alignment eltype; got $(Int)."),
                    (
                        (error("must not iterate") for _ = 1:1),
                        "Iterator must declare an alignment eltype.",
                    ),
                )
                for (input, message) in invalid_inputs
                    path = joinpath(dir, "invalid")
                    @test_throws ArgumentError(message) write_file(path, input, format)
                    @test !ispath(path)
                    write(path, "keep this file")
                    @test_throws ArgumentError(message) write_file(path, input, format)
                    @test read(path, String) == "keep this file"
                    rm(path)
                    io = IOBuffer()
                    @test_throws ArgumentError(message) print_file(io, input, format)
                    @test position(io) == 0
                end
            end
        end
    end

    @testset "Pfam example preserves mappings" begin
        fixture = joinpath(DATA, "PF09645_full.stockholm")
        expected = read_file(fixture, Stockholm)
        mktempdir() do dir
            source = joinpath(dir, "Pfam-A.full.gz")
            output = joinpath(dir, "Pfam-A.aligned.stockholm.gz")
            _write_msa_fixture(source, repeat(read(fixture, String) * "\n", 2))
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
