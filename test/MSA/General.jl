@testset "Test using MSA files" begin

    msa_types = (
        Matrix{Residue},
        NamedResidueMatrix{Array{Residue,2}},
        MultipleSequenceAlignment,
        AnnotatedMultipleSequenceAlignment,
    )

    pf09645_sto = joinpath(DATA, "PF09645_full.stockholm")
    gaoetal2011 = joinpath(DATA, "Gaoetal2011.fasta")

    gaoetal_msas = [read_file(gaoetal2011, FASTA, T) for T in msa_types]
    pfam_msas = [read_file(pf09645_sto, Stockholm, T) for T in msa_types]

    @testset "getindex" begin

        residues = permutedims(
            hcat(
                res"DAWAEE",
                res"DAWAEF",
                res"DAWAED",
                res"DAYCMD",
                res"DAYCMT",
                res"DAYCMT",
            ),
            [2, 1],
        )

        for index in [
            (2:4, 2:4),
            (:, [1, 2, 3, 4, 5, 6] .< 4),
            ([1, 2, 3, 4, 5, 6] .< 4, :),
            (:, :),
            ([1, 2, 3, 4, 5, 6] .< 4, [1, 2, 3, 4, 5, 6] .< 4),
        ]

            for msa in gaoetal_msas
                selection = msa[index...]
                @test selection == residues[index...]
                @test selection isa typeof(msa)
            end
        end

        for index in [2, (2, :), (:, 2), (2, 2), (2, [3, 4, 5])]

            for msa in gaoetal_msas
                @test msa[index...] == residues[index...]
            end
        end
    end

    @testset "setindex! and copy" begin

        for msa in gaoetal_msas
            copy_msa = copy(msa)
            deepcopy_msa = deepcopy(msa)

            for (index, value) in [((1), Residue('H')), ((:, 1), res"HHHHHH")]

                deepcopy_msa[index...] = value
                @test deepcopy_msa[index...] == value
                @test msa[index...] != value

                copy_msa[index...] = value
                @test copy_msa[index...] == value
                @test msa[index...] != value
            end

            for (index, value) in [(1, Residue('H')), (4, Residue('X')), (:, res"HHHHHH")]

                seq = copy(getsequence(msa, 4))
                seq[1, index] = value
                @test seq[1, index] == value
                @test msa[4, index] != value # since seq is a copy
            end
        end
    end

    @testset "Size" begin

        for aln in gaoetal_msas
            @test size(aln) == (6, 6)
            @test length(aln) == 36
            @test ncolumns(aln) == 6
            @test nsequences(aln) == 6
            @test ncolumns(getsequence(aln, 4)) == 6
        end

        for aln in pfam_msas
            @test size(aln) == (4, 110)
            @test length(aln) == 440
            @test ncolumns(aln) == 110
            @test nsequences(aln) == 4
        end
    end

    @testset "AnnotatedAlignedSequence and AlignedSequence" begin

        @testset "Creation" begin
            seq_types = (
                Matrix{Residue},
                NamedResidueMatrix{Array{Residue,2}},
                AlignedSequence,
                AnnotatedAlignedSequence,
            )

            for i in eachindex(msa_types)
                M = msa_types[i]
                S = seq_types[i]
                msa = pfam_msas[i]

                if M != Matrix{Residue}
                    for id in [
                        "C3N734_SULIY/1-95",
                        "H2C869_9CREN/7-104",
                        "Y070_ATV/2-70",
                        "F112_SSV1/3-112",
                    ]
                        annseq = getsequence(msa, id)
                        @test msa[id, :] == vec(annseq) # Sequences are matrices
                        @test isa(annseq, S)
                    end
                end

                for seq = 1:4
                    annseq = getsequence(msa, seq)
                    @test msa[seq, :] == vec(annseq) # Sequences are matrices
                    @test isa(annseq, S)
                end
            end
        end

        @testset "Annotations" begin
            msa = read_file(pf09645_sto, Stockholm)

            @test getannotcolumn(msa, "SS_cons") ==
                  getannotcolumn(getsequence(msa, 4), "SS_cons")
            @test getannotfile(msa) == getannotfile(getsequence(msa, 4))

            # The sequence name is only needed when working with MSA objects.
            @test getannotresidue(msa, "F112_SSV1/3-112", "SS") ==
                  getannotresidue(getsequence(msa, 4), "SS")
            @test getannotsequence(msa, "F112_SSV1/3-112", "DR") ==
                  getannotsequence(getsequence(msa, 4), "DR")
        end
    end

    @testset "Print" begin
        pfam = read_file(pf09645_sto, Stockholm)
        gao = read_file(gaoetal2011, FASTA)

        @test stringsequence(gao, 4) == "DAYCMD"
        @test stringsequence(pfam, 1) == stringsequence(pfam, "C3N734_SULIY/1-95")

        for T in (Stockholm, FASTA, Raw)
            buffer = IOBuffer()
            print_file(buffer, pfam, T)
            @test parse_file(String(take!(buffer)), T) == pfam
            print_file(buffer, gao, T)
            @test parse_file(String(take!(buffer)), T) == gao
        end
    end
end

@testset "Base.show MIME text/plain for MSA Containers" begin
    using MIToS.MSA
    using MIToS.Utils # For NamedArray, All
    using NamedArrays
    using OrderedCollections # For OrderedDict in Annotations

    # Helper to create a simple NamedArray for sequences
    function create_simple_named_array(nseq, ncol)
        data = rand(Residue, nseq, ncol)
        seqnames = ["Seq" * string(i) for i in 1:nseq]
        colnames = ["Col" * string(j) for j in 1:ncol]
        return NamedArray(data, (OrderedDict(k => i for (i,k) in enumerate(seqnames)), 
                                 OrderedDict(k => j for (j,k) in enumerate(colnames))), ("Sequences", "Columns"))
    end

    @testset "MultipleSequenceAlignment show" begin
        na = create_simple_named_array(2, 3)
        msa = MultipleSequenceAlignment(na)
        
        # Using sprint with MIME type
        str_output = sprint((io, x) -> show(io, MIME"text/plain"(), x), msa)
        
        @test startswith(str_output, "MultipleSequenceAlignment : ")
        @test occursin("2x3 Named Matrix{Residue", str_output) # Part of NamedArray show
        @test occursin("Seq1", str_output)
        @test occursin("Col1", str_output)
    end

    @testset "AnnotatedMultipleSequenceAlignment show" begin
        na = create_simple_named_array(2, 3)
        msa_part = MultipleSequenceAlignment(na)
        
        file_annot = OrderedDict("GF_ID" => ["GF Value"])
        seq_annot = OrderedDict("Seq1" => OrderedDict("GS_ID" => ["GS Seq1 Value"]))
        col_annot = OrderedDict("Col1" => OrderedDict("GC_ID" => ["GC Col1 Value"]))
        res_annot = OrderedDict(("Seq1", "Col1") => OrderedDict("GR_ID" => ["GR S1C1 Value"]))
        
        annotations = MSAAnnotations(
            file_annot,
            seq_annot,
            col_annot,
            res_annot
        )
        annot_msa = AnnotatedMultipleSequenceAlignment(msa_part, annotations)

        str_output = sprint((io, x) -> show(io, MIME"text/plain"(), x), annot_msa)
        
        # Count total distinct annotation keys (GF_ID, GS_ID, GC_ID, GR_ID)
        # This is a bit tricky as the "with X annotations" refers to annotation groups/categories.
        # Based on `_show_annotations_MIToS`, it counts the number of keys in each dictionary (file, seq, col, res)
        # and sums them if they are non-empty.
        # Let's assume "with X annotations" refers to the number of distinct annotation types present.
        # The current implementation counts the number of annotation *categories* that have entries.
        # So if file, seq, col, res all have entries, it's "with 4 annotations"
        # If only file and seq have entries, it's "with 2 annotations"
        # Given the example from problem description: "AnnotatedMultipleSequenceAlignment with Y annotations : "
        # This means Y is the number of annotation *groups* that are non-empty.
        # Here, all 4 (file, sequence, column, residue) have annotations.
        @test startswith(str_output, "AnnotatedMultipleSequenceAlignment with 4 annotations : ")
        @test occursin("2x3 Named Matrix{Residue", str_output)
        @test occursin("Seq1", str_output)
        @test occursin("Col1", str_output)

        # Test with fewer annotations
        annotations_minimal = MSAAnnotations(file_annot, OrderedDict{String,OrderedDict{String,Vector{String}}}()) # only file annotations
        annot_msa_minimal = AnnotatedMultipleSequenceAlignment(msa_part, annotations_minimal)
        str_output_minimal = sprint((io, x) -> show(io, MIME"text/plain"(), x), annot_msa_minimal)
        @test startswith(str_output_minimal, "AnnotatedMultipleSequenceAlignment with 1 annotations : ")
    end

    @testset "AlignedSequence show" begin
        na = create_simple_named_array(1, 5) # AlignedSequence is effectively 1xN
        aligned_seq_data = MultipleSequenceAlignment(na)[1,:] # Get an AlignedSequence
        
        str_output = sprint((io, x) -> show(io, MIME"text/plain"(), x), aligned_seq_data)
        
        @test startswith(str_output, "AlignedSequence : ")
        @test occursin("1x5 Named Matrix{Residue", str_output) # From NamedArray
        @test occursin("Seq1", str_output) # Sequence name
        @test occursin("Col1", str_output) # Column name
    end

    @testset "AnnotatedAlignedSequence show" begin
        na = create_simple_named_array(1, 5)
        msa_part = MultipleSequenceAlignment(na)
        
        file_annot = OrderedDict("GF_ID" => ["GF Value"]) # These won't show directly for sequence
        # Sequence specific annotations for "Seq1"
        seq_annot_for_seq1 = OrderedDict("Seq1" => OrderedDict("GS_ID" => ["GS Seq1 Value"]))
        # Column specific annotations for "Col1"
        col_annot_for_col1 = OrderedDict("Col1" => OrderedDict("GC_ID" => ["GC Col1 Value"]))
        # Residue specific for ("Seq1", "Col1")
        res_annot_for_s1c1 = OrderedDict(("Seq1", "Col1") => OrderedDict("GR_ID" => ["GR S1C1 Value"]))

        annotations = MSAAnnotations(file_annot, seq_annot_for_seq1, col_annot_for_col1, res_annot_for_s1c1)
        annot_msa = AnnotatedMultipleSequenceAlignment(msa_part, annotations)
        
        # Get the AnnotatedAlignedSequence
        annot_aligned_seq = getsequence(annot_msa, "Seq1") # Or annot_msa[1,:]

        str_output = sprint((io, x) -> show(io, MIME"text/plain"(), x), annot_aligned_seq)
        
        # AnnotatedAlignedSequence show method refers to sequence, column and residue annotations for that sequence.
        # It sums the count of (sequence specific), (column specific for that seq), (residue specific for that seq)
        # Here: GS_ID (1), GC_ID (1 for Col1), GR_ID (1 for S1C1) -> 3 annotation types for this sequence
        # The actual number of annotation *groups* associated with this sequence.
        # The `_show_annotations_MIToS` for sequence counts non-empty sequence, column, and residue annotations.
        # For this sequence "Seq1":
        # - It has sequence annotations (GS_ID).
        # - It has column annotations (GC_ID for Col1).
        # - It has residue annotations (GR_ID for (Seq1,Col1)).
        # So, 3 groups of annotations are relevant.
        @test startswith(str_output, "AnnotatedAlignedSequence with 3 annotations : ")
        @test occursin("1x5 Named Matrix{Residue", str_output)
        @test occursin("Seq1", str_output)
        @test occursin("Col1", str_output)

        # Test with only sequence annotation for this sequence
        annotations_seq_only = MSAAnnotations(OrderedDict{String,Vector{String}}(), seq_annot_for_seq1)
        annot_msa_seq_only = AnnotatedMultipleSequenceAlignment(msa_part, annotations_seq_only)
        annot_aligned_seq_minimal = getsequence(annot_msa_seq_only, "Seq1")
        str_output_minimal = sprint((io, x) -> show(io, MIME"text/plain"(), x), annot_aligned_seq_minimal)
        @test startswith(str_output_minimal, "AnnotatedAlignedSequence with 1 annotations : ")
    end
    
    @testset "AnnotatedSequence show" begin
        # AnnotatedSequence is a bit different. It's a general container for a sequence and its annotations.
        # It does not necessarily come from an MSA and might not have column names in the same way.
        # Let's create one directly.
        seq_data = res"ACGT"
        seq_name = "MySeq"
        # Annotations for this sequence (e.g., features, comments)
        # These are not per-residue or per-column in the MSA sense, but rather for the sequence itself.
        # The structure of annotations for AnnotatedSequence might be simpler, e.g. a Dict or NamedTuple.
        # Looking at `src/MSA/AnnotatedSequence.jl`, it stores `sequence::AbstractVector{Residue}`
        # and `annotations::Annotations`. The `Annotations` here is the generic one from `src/MSA/Annotations.jl`
        # which can store GF, GS, GC, GR. For an isolated sequence, GS would be most relevant.
        
        seq_annotations_dict = OrderedDict{String,Vector{String}}("Feature" => ["Protein kinase domain"], "Source" => ["Manual annotation"])
        # Wrap this in the main Annotations struct, possibly as file annotations if no other context
        # Or, more likely, these are sequence-level annotations if AnnotatedSequence is to be part of an MSA.
        # If it's an isolated sequence, these are just metadata.
        # The `show` method for `AnnotatedSequence` currently prints:
        # `AnnotatedSequence with X annotations : ID : Sequence`
        # X is the number of sequence-level annotations for its ID.
        
        # Create Annotations that would apply to this sequence if it had an ID "MySeq"
        gs_annotations = MSAAnnotations( # File annots
                                        OrderedDict{String,Vector{String}}(), 
                                        # Sequence annots
                                        OrderedDict("MySeq" => seq_annotations_dict),
                                        # Column annots (empty for isolated seq)
                                        OrderedDict{String,OrderedDict{String,Vector{String}}}()) 
                                        # Residue annots (empty for isolated seq)

        # Create an AnnotatedSequence. It needs an ID.
        # The constructor is `AnnotatedSequence(id::String, sequence::AbstractVector{Residue}, annotations::Annotations)`
        # The `annotations` here is the global MSAAnnotations object.
        annot_seq = AnnotatedSequence("MySeq", seq_data, gs_annotations)

        str_output = sprint((io, x) -> show(io, MIME"text/plain"(), x), annot_seq)

        # "MySeq" has 2 sequence annotations ("Feature", "Source")
        @test startswith(str_output, "AnnotatedSequence with 2 annotations : MySeq : ")
        @test occursin("A C G T", str_output) # Default printing of sequence vector

        # Test with no specific annotations for this ID
        empty_annotations = MSAAnnotations()
        annot_seq_no_annot = AnnotatedSequence("NoAnnotSeq", seq_data, empty_annotations)
        str_output_no_annot = sprint((io, x) -> show(io, MIME"text/plain"(), x), annot_seq_no_annot)
        @test startswith(str_output_no_annot, "AnnotatedSequence with 0 annotations : NoAnnotSeq : ")
        @test occursin("A C G T", str_output_no_annot)
    end
end
