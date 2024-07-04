@testset "AnnotatedSequence" begin

    seq = AnnotatedSequence("id", res"ARNDCQEGHILKMFPSTWYV", Annotations())
    seq_from_string = AnnotatedSequence("id", "arn-.DCQEGHILKMFPSTWYV", Annotations())
    seq_no_id = AnnotatedSequence(res"ARNDCQEGHILKMFPSTWYV", Annotations())
    seq_from_string_no_id = AnnotatedSequence("arn-.DCQEGHILKMFPSTWYV", Annotations())

    sequences = [seq, seq_from_string, seq_no_id, seq_from_string_no_id]

    @testset "Creation" begin

        for sequence in sequences
            @test dimnames(sequence) == ["Seq", "Pos"]
            @test namedmatrix(sequence) == namedmatrix(seq) # as the ids are different
        end
    end

    @testset "Interfaces : AbstractArray" begin

        @test IndexStyle(AnnotatedSequence) == IndexLinear()

        @test size(seq) == (1, 20) # A sequence is stored as a 1xN matrix
        @test length(seq) == 20
        @test getindex(seq, 1) == Residue('A')
        @test getindex(seq, 1:3) == res"ARN"
        @test join(Char(res) for res in seq) == "ARNDCQEGHILKMFPSTWYV" # Iteration
    end

    @testset "Sequence ID" begin

        @test sequence_id(seq) == "id"
        @test sequence_id(seq_from_string) == "id"
        @test sequence_id(seq_no_id) == ""
        @test sequence_id(seq_from_string_no_id) == ""
    end

    @testset "Equality" begin

        # ==
        @test seq == seq_from_string # same id and sequence
        @test seq_no_id == seq_from_string_no_id
        @test seq !== seq_no_id # different id, same sequence
        @test seq_from_string !== seq_from_string_no_id

        # hash
        @test hash(seq) == hash(seq_from_string)
        @test hash(seq_no_id) == hash(seq_from_string_no_id)
        @test hash(seq) !== hash(seq_no_id)
        @test hash(seq_from_string) !== hash(seq_from_string_no_id)

        # isequal
        @test isequal(seq, seq_from_string)
        @test isequal(seq_no_id, seq_from_string_no_id)
        @test !isequal(seq, seq_no_id)
        @test !isequal(seq_from_string, seq_from_string_no_id)
    end

    @testset "String" begin

        @test stringsequence(seq) == "ARNDCQEGHILKMFPSTWYV"
        @test join(seq) == "ARNDCQEGHILKMFPSTWYV"
    end

    @testset "From file" begin
        pir = read_file(joinpath(DATA, "emboss.pir"), PIRSequences)

        @test isa(pir, Vector{AnnotatedSequence})
        @test length(pir) == 1 # there is only one sequence in the file

        seq = pir[1]
        @test join(seq[1:12]) == "GDVEGKGIFTMC"
        @test sequence_id(seq) == "AZBR"

        @test getannotsequence(seq, "Type") == "P1"
        @test getannotsequence(seq, "Title") ==
              "finger protein zfpA - turnip fern chloroplast"

        # This should show a warning, as the sequence name is given as a feature
        @test getannotsequence(seq, "AZBR", "Type") == "Type" # returns the default value
        @test getannotsequence(seq, "AZBR", "Title") == "Title"

        @testset "ThorAxe's transcripts.pir" begin
            # This file contains the first four sequences from the ThorAxe output for POLR3B
            # The annotations correspond to s-exon symbols. There is one symbol per residue.
            # The list of symbols is at the end of the sequence identifier.
            transcripts = read_file(joinpath(DATA, "transcripts.pir"), PIRSequences)

            @test isa(transcripts, Vector{AnnotatedSequence})
            @test length(transcripts) == 4

            for seq in transcripts
                @test !isempty(annotations(seq))
                @test length(seq) == length(getannotsequence(seq, "Title"))
                @test split(sequence_id(seq))[end] ==
                      join(unique(getannotsequence(seq, "Title")))
            end
        end

        @testset "UniProt's FASTA (canonical & isoform)" begin
            # This file contains the canonical and isoform sequences for the human protein
            # POLR3B (Q9NW08) in FASTA format as downloaded from UniProt.
            fasta = read_file(joinpath(DATA, "Q9NW08.fasta"), FASTASequences)

            @test isa(fasta, Vector{AnnotatedSequence})
            @test length(fasta) == 2

            Ns = [1_133, 1_075]
            for (i, seq) in enumerate(fasta)
                @test startswith(sequence_id(seq), "sp|Q9NW08")
                @test isempty(getannotsequence(seq))
                @test endswith(join(seq), "SKYNE")
                @test length(seq) == Ns[i]
            end
        end
    end

    @testset "From an MSA" begin
        raw = read_file(joinpath(DATA, "gaps.txt"), RawSequences)
        fasta = read_file(joinpath(DATA, "ids.fasta"), FASTASequences) # same sequence, different ids

        @test isa(raw, Vector{AnnotatedSequence})
        @test isa(fasta, Vector{AnnotatedSequence})

        @test length(raw) == 10
        @test length(fasta) == 4

        for (i, seq) in enumerate(raw) # from gaps.txt
            n_res = 10 - (i - 1) # the first sequence has 0 gaps, and the last one has 9 gaps
            @test size(seq) == (1, n_res)
            @test length(seq) == n_res
            @test sequence_id(seq) == string(i)
            @test isempty(annotations(seq)) # no annotations
        end

        for i = 2:4
            @test join(fasta[i]) == join(fasta[1]) # all sequences have the same residues
        end
        # but different ids
        @test sequence_id.(fasta) == [
            "SEQUENCE_1",
            " SEQUENCE_2",
            "MCHU - Calmodulin - Human, rabbit, bovine, rat, and chicken",
            "gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]",
        ]
    end

    @testset "Getter Functions" begin
        annot = Annotations(
            OrderedDict("FileFeature" => "FileAnnotation"),
            Dict(("SeqID", "Type") => "P1", ("SeqID", "Title") => "Protein Title"),
            Dict("ColFeature" => "HHHHHH"),
            Dict(("SeqID", "ResFeature") => "001100"),
        )

        seq = AnnotatedSequence("SeqID", "ARNDCQ", annot)

        @test getannotfile(seq) == OrderedDict("FileFeature" => "FileAnnotation")
        @test getannotcolumn(seq) == Dict("ColFeature" => "HHHHHH")
        @test getannotsequence(seq, "Type") == "P1"
        @test getannotsequence(seq, "Title") == "Protein Title"
        @test getannotresidue(seq, "ResFeature") == "001100"
    end

    @testset "Setter Functions" begin
        new_seq = deepcopy(seq)
        setannotfile!(new_seq, "NewFileFeature", "NewFileAnnotation")
        @test getannotfile(new_seq, "NewFileFeature") == "NewFileAnnotation"

        setannotcolumn!(new_seq, "NewColFeature", "::..  ")
        @test getannotcolumn(new_seq, "NewColFeature") == "::..  "

        setannotsequence!(new_seq, "NewType", "NewTypeAnnotation")
        @test getannotsequence(new_seq, "NewType") == "NewTypeAnnotation"

        setannotresidue!(new_seq, "NewResFeature", "ARnd..")
        @test getannotresidue(new_seq, "NewResFeature") == "ARnd.."
    end

    @testset "IO" begin
        @testset "FASTASequences" begin
            in_fasta = read_file(joinpath(DATA, "Q9NW08.fasta"), FASTASequences)
            io = IOBuffer()
            print_file(io, in_fasta, FASTASequences)
            seekstart(io)
            out_fasta = parse_file(io, FASTASequences)
            @test in_fasta == out_fasta
        end

        @testset "PIRSequences" begin
            in_pir = read_file(joinpath(DATA, "transcripts.pir"), PIRSequences)
            io = IOBuffer()
            print_file(io, in_pir, PIRSequences)
            seekstart(io)
            out_pir = parse_file(io, PIRSequences)
            @test in_pir == out_pir
        end

        @testset "RawSequences" begin
            in_raw = read_file(joinpath(DATA, "raw_sequences.txt"), RawSequences)
            io = IOBuffer()
            print_file(io, in_raw, RawSequences)
            seekstart(io)
            out_raw = parse_file(io, RawSequences)
            @test in_raw == out_raw
        end
    end
end
