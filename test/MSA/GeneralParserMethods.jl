@testset "GeneralParserMethods" begin

    @testset "Test input lengths" begin
        # NOTE: _pre_read... functions should call _check_seq_len
        # if _convert_to_matrix_residues is used

        ids = ["A", "B", "C"]
        seqs_equal = ["MSRSKRDNEIGDSTF", "MSRSKRDNNFGDSTF", "MSRSFYSVEIGDSTF"]
        seqs_diffs = ["MSRSKRDNEIGDSTF", "MSRSKRDNNFGDSTF", "MSRSFYSVEIGD"]
        seqs_less = ["MSRSKRDNEIGDSTF", "MSRSKRDNNFGDSTF"]

        @test MSA._check_seq_len(ids, seqs_equal) === nothing
        @test_throws ErrorException MSA._check_seq_len(ids, seqs_diffs)
        @test_throws ErrorException MSA._check_seq_len(ids, seqs_less)
    end

    @testset "To parse MSA & mapping" begin

        sequences = ["ADEIMSY", "RCGLFTV", "NQHKPW-"]
        matrixres = reshape(reinterpret(Residue, collect(1:21)), (3, 7))
        msa, map = MSA._to_msa_mapping(sequences)

        @test getarray(msa) == matrixres
        @test map == ["1,2,3,4,5,6,7", "1,2,3,4,5,6,7", "1,2,3,4,5,6,"]
        @test sequencenames(msa) == ["1", "2", "3"]
        # MSA constructors adds dimension names
        # @test dimnames(msa) == [:Seq, :Col]

        msa, map = MSA._to_msa_mapping(sequences, ["a/11-17", "b/11-17", "c/11-16"])

        @test getarray(msa) == matrixres
        @test map == ["11,12,13,14,15,16,17", "11,12,13,14,15,16,17", "11,12,13,14,15,16,"]
        @test sequencenames(msa) == ["a/11-17", "b/11-17", "c/11-16"]
        # MSA constructors adds dimension names
        # @test dimnames(msa) == [:Seq, :Col]

        @test_throws ErrorException MSA._to_msa_mapping(
            sequences,
            ["a/1-7", "b/1-7", "c/1-7"],
        )
        @test_throws ErrorException MSA._to_msa_mapping(["AD", "EI", "MSY"])
    end

    @testset "Delete full gap columns" begin

        M = reshape(reinterpret(Residue, collect(1:21)), (3, 7))
        M[:, [2, 4, 6]] .= GAP
        msa = MultipleSequenceAlignment(M)
        named = namedmatrix(msa)
        annotated_msa = AnnotatedMultipleSequenceAlignment(M)

        d_M = deletefullgapcolumns(M)
        d_msa = deletefullgapcolumns(msa)
        d_named = deletefullgapcolumns(named)
        d_annotated_msa = deletefullgapcolumns(annotated_msa)

        @test size(d_M) == (3, 4)

        for object in [d_msa, d_named, d_annotated_msa]
            @test size(object) == (3, 4)
            @test getcolumnmapping(object) == [1, 3, 5, 7]
        end

        d_msa = deletefullgapcolumns!(msa)
        d_annotated_msa = deletefullgapcolumns!(annotated_msa)

        @test d_msa == msa
        @test d_annotated_msa == annotated_msa
    end


    @testset "Disambiguate Sequences (online) Tests" begin
        # helper that returns both the disambiguator and the list of new IDs
        function run_online(ids)
            disambiguator = MSA.OnlineSequenceNameDisambiguator()
            new_ids = [MSA._disambiguate_seqname!(disambiguator, id) for id in ids]
            return disambiguator, new_ids
        end

        run_bulk(ids) = MSA._disambiguate_seqnames!(copy(ids), Annotations()) # new_ids, annot

        ## Test 1: single raw_id repeated
        ids = ["seq", "seq", "seq"]
        disambiguator, new_ids = run_online(ids)
        @test new_ids == ["seq", "seq(1)", "seq(2)"]
        @test disambiguator.new2old ==
              OrderedDict("seq" => "seq", "seq(1)" => "seq", "seq(2)" => "seq")
        new_ids, annot = run_bulk(ids)
        @test new_ids == ["seq", "seq(1)", "seq(2)"]
        @test getannotsequence(annot, "seq", "OriginalSeqName", "") == ""
        @test getannotsequence(annot, "seq(1)", "OriginalSeqName", "") == "seq"
        @test getannotsequence(annot, "seq(2)", "OriginalSeqName", "") == "seq"

        ## Test 2: two raw_ids interleaved
        ids = ["a", "a", "b", "b", "b"]
        disambiguator, new_ids = run_online(ids)
        @test new_ids == ["a", "a(1)", "b", "b(1)", "b(2)"]
        @test disambiguator.new2old == OrderedDict(
            "a" => "a",
            "a(1)" => "a",
            "b" => "b",
            "b(1)" => "b",
            "b(2)" => "b",
        )
        new_ids, annot = run_bulk(ids)
        @test new_ids == ["a", "a(1)", "b", "b(1)", "b(2)"]
        @test getannotsequence(annot, "a", "OriginalSeqName", "") == ""
        @test getannotsequence(annot, "a(1)", "OriginalSeqName", "") == "a"
        @test getannotsequence(annot, "b", "OriginalSeqName", "") == ""
        @test getannotsequence(annot, "b(1)", "OriginalSeqName", "") == "b"
        @test getannotsequence(annot, "b(2)", "OriginalSeqName", "") == "b"

        ## Test 3: names that already contain suffixes
        ids = ["item", "item(1)", "item(1)", "item(2)"]
        disambiguator, new_ids = run_online(ids)
        @test new_ids == ["item", "item(1)", "item(1)(1)", "item(2)"]
        @test disambiguator.new2old == OrderedDict(
            "item" => "item",
            "item(1)" => "item(1)",
            "item(1)(1)" => "item(1)",
            "item(2)" => "item(2)",
        )
        new_ids, annot = run_bulk(ids)
        @test new_ids == ["item", "item(1)", "item(1)(1)", "item(2)"]
        @test getannotsequence(annot, "item", "OriginalSeqName", "") == ""
        @test getannotsequence(annot, "item(1)", "OriginalSeqName", "") == ""
        @test getannotsequence(annot, "item(1)(1)", "OriginalSeqName", "") == "item(1)"
        @test getannotsequence(annot, "item(2)", "OriginalSeqName", "") == ""

        ## Test 4: nested suffix collision
        ids = ["x", "x", "x", "x(1)", "x(1)"]
        disambiguator, new_ids = run_online(ids)
        @test new_ids == ["x", "x(1)", "x(2)", "x(1)(1)", "x(1)(2)"]
        @test disambiguator.new2old == OrderedDict(
            "x" => "x",
            "x(1)" => "x",
            "x(2)" => "x",
            "x(1)(1)" => "x(1)",
            "x(1)(2)" => "x(1)",
        )
        new_ids, annot = run_bulk(ids)
        @test new_ids == ["x", "x(1)", "x(2)", "x(1)(1)", "x(1)(2)"]
        @test getannotsequence(annot, "x", "OriginalSeqName", "") == ""
        @test getannotsequence(annot, "x(1)", "OriginalSeqName", "") == "x"
        @test getannotsequence(annot, "x(2)", "OriginalSeqName", "") == "x"
        @test getannotsequence(annot, "x(1)(1)", "OriginalSeqName", "") == "x(1)"
        @test getannotsequence(annot, "x(1)(2)", "OriginalSeqName", "") == "x(1)"
    end

    @testset "Mapping warnings" begin
        fasta = joinpath(DATA, "simple.fasta")
        msa = read_file(fasta, FASTA, generatemapping = true)
        mktemp() do tmp, _
            write_file(tmp, msa, Stockholm)
            @test_logs (:warn, r"sequence mappings annotations") (
                :warn,
                r"column annotations",
            ) match_mode = :all begin
                read_file(tmp, Stockholm, generatemapping = true)
            end
        end
    end
end
