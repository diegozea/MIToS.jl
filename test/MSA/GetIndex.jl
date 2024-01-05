@testset "getindex" begin
    simple = joinpath(DATA, "simple.fasta")

    @testset "get index" begin
        msa = read(simple, FASTA, generatemapping=true)
        first_seq = getsequence(msa, 1)
        matrix_msa = convert(Matrix{Residue}, msa)

        @testset "sequence_index" begin
            @test sequence_index(msa, "ONE") == 1
            @test sequence_index(msa, "TWO") == 2
            @test_throws KeyError sequence_index(msa, "THREE")  # non-existent sequence

            # Matrix{Residue}
            @test_throws ErrorException sequence_index(matrix_msa, "ONE")

            # AlignedSequence
            @test first_seq isa AbstractAlignedSequence
            @test sequence_index(first_seq, "ONE") == 1
            @test_throws KeyError sequence_index(first_seq, "TWO")  # non-existent sequence

            # Test returning the same index when an integer is passed
            @test sequence_index(msa, 1) == 1
            @test sequence_index(msa, 2) == 2
        end

        @testset "column_index" begin
            @test column_index(msa, "1") == 1
            @test column_index(msa, "2") == 2
            @test_throws KeyError column_index(msa, "3")  # non-existent column

            # Matrix{Residue}
            @test_throws ErrorException column_index(matrix_msa, "1")

            # AlignedSequence
            @test first_seq isa AbstractAlignedSequence
            @test column_index(first_seq, "1") == 1
            @test column_index(first_seq, "2") == 2
            @test_throws KeyError column_index(first_seq, "3")  # non-existent column

            # Test returning the same index when an integer is passed
            @test column_index(msa, 1) == 1
            @test column_index(msa, 2) == 2
        end
    end

    @testset "MSA" begin
        msa = read(simple, FASTA, generatemapping=true)

        matrix = Residue['A'  'R'
                         'R'  'A']
        
        @test msa == matrix
        @test getcolumnmapping(msa) == [1, 2]
        @test getsequencemapping(msa, "ONE") == [1, 2]
        @test getsequencemapping(msa, "TWO") == [1, 2]

        for reverse_col_selector in [[2, 1], 2:-1:1, String["2", "1"]]
            ref = Residue['R'  'A'; 'A'  'R']

            reversed_cols = msa[:, reverse_col_selector]
            @test reversed_cols == ref
            @test getcolumnmapping(reversed_cols) == [2, 1]
            @test getsequencemapping(reversed_cols, "ONE") == [2, 1]
            @test getsequencemapping(reversed_cols, "TWO") == [2, 1]
            @test sequencenames(reversed_cols) == ["ONE", "TWO"]
            @test columnnames(reversed_cols) == ["2", "1"]
            annot_values = Set(values(getannotfile(reversed_cols)))
            @test any(occursin("filtercolumns! : 2 columns have been selected.", val) for val in annot_values)
            @test any(occursin("filtercolumns! : column order has changed!", val) for val in annot_values)

            @test MultipleSequenceAlignment(msa)[:, reverse_col_selector] == ref
        end

        for reverse_seq_selector in [[2, 1], 2:-1:1, String["TWO", "ONE"]]
            ref = Residue['R'  'A'; 'A'  'R']

            reversed_seqs = msa[reverse_seq_selector, :]
            @test reversed_seqs == ref
            @test getcolumnmapping(reversed_seqs) == [1, 2]
            @test getsequencemapping(reversed_seqs, "ONE") == [1, 2]
            @test getsequencemapping(reversed_seqs, "TWO") == [1, 2]
            @test sequencenames(reversed_seqs) == ["TWO", "ONE"]
            @test columnnames(reversed_seqs) == ["1", "2"]

            @test MultipleSequenceAlignment(msa)[reverse_seq_selector, :] == ref
        end

        ref_single_col = matrix[:, [2]]

        for single_col_selector in [[2], [false, true], String["2"]]
            single_col = msa[:, single_col_selector]
            @test single_col == ref_single_col
            @test getcolumnmapping(single_col) == [2]
            @test getsequencemapping(single_col, "ONE") == [2]
            @test getsequencemapping(single_col, "TWO") == [2]
            @test sequencenames(single_col) == ["ONE", "TWO"]
            @test columnnames(single_col) == ["2"]
            annot_values = Set(values(getannotfile(single_col)))
            @test any(occursin("filtercolumns! : 1 column has been", val) for val in annot_values)

            @test MultipleSequenceAlignment(msa)[:, single_col_selector] == ref_single_col
        end

        ref_single_seq = matrix[[2], :]

        for single_seq_selector in [[2], [false, true], String["TWO"]]
            single_seq = msa[single_seq_selector, :]
            @test single_seq == ref_single_seq
            @test getcolumnmapping(single_seq) == [1, 2]
            @test_throws KeyError getsequencemapping(single_seq, "ONE")
            @test getsequencemapping(single_seq, "TWO") == [1, 2]
            @test sequencenames(single_seq) == ["TWO"]
            @test columnnames(single_seq) == ["1", "2"]
            annot_values = Set(values(getannotfile(single_seq)))
            @test any(occursin("filtersequences! : 1 sequence has been", val) for val in annot_values)

            @test MultipleSequenceAlignment(msa)[single_seq_selector, :] == ref_single_seq
        end

        ref_single_res = matrix[[2], [2]]

        for single_seq_selector in [[2], [false, true], String["TWO"]]
            for single_col_selector in [[2], [false, true], String["2"]]
                single_res = msa[single_seq_selector, single_col_selector]
                @test single_res == ref_single_res
                @test getcolumnmapping(single_res) == [2]
                @test_throws KeyError getsequencemapping(single_res, "ONE")
                @test getsequencemapping(single_res, "TWO") == [2]
                @test sequencenames(single_res) == ["TWO"]
                @test columnnames(single_res) == ["2"]
                annot_values = Set(values(getannotfile(single_res)))
                @test any(occursin("filtersequences! : 1 sequence has been", val) for val in annot_values)
                @test any(occursin("filtercolumns! : 1 column has been", val) for val in annot_values)

                @test MultipleSequenceAlignment(msa)[single_seq_selector, single_col_selector] == ref_single_res
            end
        end
    end
end