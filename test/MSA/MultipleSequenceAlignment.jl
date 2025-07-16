@testset "MultipleSequenceAlignment" begin

    @testset "Type Hierarchy" begin

        @test AbstractAlignedObject <: AbstractMatrix{Residue}
        @test AbstractMultipleSequenceAlignment <: AbstractAlignedObject
        @test AbstractAlignedSequence <: AbstractAlignedObject
        @test MultipleSequenceAlignment <: AbstractMultipleSequenceAlignment
        @test AnnotatedMultipleSequenceAlignment <: AbstractMultipleSequenceAlignment
        @test AlignedSequence <: AbstractAlignedSequence
        @test AnnotatedAlignedSequence <: AbstractAlignedSequence

        @test AnnotatedMultipleSequenceAlignment <: AnnotatedAlignedObject
        @test AnnotatedAlignedSequence <: AnnotatedAlignedObject

        @test MultipleSequenceAlignment <: UnannotatedAlignedObject
        @test AlignedSequence <: UnannotatedAlignedObject
    end

    @testset "Creation" begin

        M = reshape(reinterpret(Residue, collect(1:21)), (3, 7))
        m = NamedArray(M)

        T = permutedims(M, [2, 1])

        S = reshape(reinterpret(Residue, collect(1:21)), (1, 21))
        s = NamedArray(S)

        @test MultipleSequenceAlignment(M) == MultipleSequenceAlignment(m)
        @test AnnotatedMultipleSequenceAlignment(M) == AnnotatedMultipleSequenceAlignment(m)

        @test AlignedSequence(S) == AlignedSequence(s)
        @test AnnotatedAlignedSequence(S) == AnnotatedAlignedSequence(s)

        @test MultipleSequenceAlignment(M) != MultipleSequenceAlignment(T)
        @test AnnotatedMultipleSequenceAlignment(M) != AnnotatedMultipleSequenceAlignment(T)

        @test_throws ArgumentError AlignedSequence(M)
        @test_throws ArgumentError AnnotatedAlignedSequence(M)
    end

    @testset "MSA & sequences" begin

        M = reshape(reinterpret(Residue, collect(1:21)), (3, 7))
        msa = MultipleSequenceAlignment(M)
        annotated_msa = AnnotatedMultipleSequenceAlignment(M)

        S = reshape(reinterpret(Residue, collect(1:21)), (1, 21))
        sequence = AlignedSequence(S)
        annotated_sequence = AnnotatedAlignedSequence(S)

        @testset "Getters" begin
            @test isa(annotations(annotated_msa), Annotations)
            @test isa(annotations(annotated_sequence), Annotations)

            # unannotated objects return empty annotations
            for unannotated_object in
                (Matrix{Residue}(M), msa, Matrix{Residue}(S), sequence)
                @test isempty(annotations(unannotated_object))
                @test isa(annotations(unannotated_object), Annotations)
            end

            for object in (msa, annotated_msa, sequence, annotated_sequence)
                @test isa(namedmatrix(object), NamedResidueMatrix{Array{Residue,2}})
            end
        end

        @testset "Dimension names" begin

            for object in (msa, annotated_msa, sequence, annotated_sequence)
                @test dimnames(namedmatrix(object)) == ["Seq", "Col"]
            end
        end

        @testset "Convert" begin

            msa2annot = AnnotatedMultipleSequenceAlignment(msa)
            annot2msa = MultipleSequenceAlignment(annotated_msa)

            seq2annot = AnnotatedAlignedSequence(sequence)
            annot2seq = AlignedSequence(annotated_sequence)

            @test isa(annotations(msa2annot), Annotations)
            @test isa(annotations(seq2annot), Annotations)

            # unannotated objects return empty annotations
            for unannotated_object in (annot2msa, annot2seq)
                @test isempty(annotations(unannotated_object))
                @test isa(annotations(unannotated_object), Annotations)
            end

            @test namedmatrix(msa2annot) == namedmatrix(msa)
            @test namedmatrix(annot2msa) == namedmatrix(annotated_msa)
            @test namedmatrix(seq2annot) == namedmatrix(sequence)
            @test namedmatrix(annot2seq) == namedmatrix(annotated_sequence)
        end

        @testset "AbstractArray Interface" begin

            @test size(msa) == (3, 7)
            @test size(annotated_msa) == (3, 7)

            @test size(sequence) == (1, 21)
            @test size(annotated_sequence) == (1, 21)

            for object in (msa, annotated_msa, sequence, annotated_sequence)
                @test length(object) == 21
            end
        end

        @testset "Indexing" begin

            for object in (msa, annotated_msa)
                for i = 1:3, j = 1:7
                    @test object[string(i), string(j)] == object[i, j]
                end
            end

            for object in (sequence, annotated_sequence)
                for j = 1:7
                    @test object["1", string(j)] == object[1, j]
                    @test object[j] == object[1, j]
                    # Special sequence indexing:
                    @test object[string(j)] == object[1, j]
                end
            end
        end

        @testset "Show" begin

            out = IOBuffer()

            show(out, MIME"text/plain"(), msa)
            str = String(take!(out))
            @test startswith(str, "MultipleSequenceAlignment : ")
            @test occursin("Seq", str)
            @test occursin("Col", str)
            @test length(split(str, '\n')) == 6

            show(out, MIME"text/plain"(), annotated_msa)
            str = String(take!(out))
            @test startswith(
                str,
                "AnnotatedMultipleSequenceAlignment with 0 annotations : ",
            )
            @test occursin("Seq", str)
            @test occursin("Col", str)
            @test length(split(str, '\n')) == 6

            show(out, MIME"text/plain"(), sequence)
            str = String(take!(out))
            @test startswith(str, "AlignedSequence : ")
            @test occursin("Seq", str)
            @test occursin("Col", str)
            @test length(split(str, '\n')) == 4

            show(out, MIME"text/plain"(), annotated_sequence)
            str = String(take!(out))
            @test startswith(str, "AnnotatedAlignedSequence with 0 annotations : ")
            @test occursin("Seq", str)
            @test occursin("Col", str)
            @test length(split(str, '\n')) == 4
        end

        @testset "Transpose, i.e. permutedims" begin

            @test size(permutedims(msa)) == (7, 3)
            @test size(permutedims(annotated_msa)) == (7, 3)

            @test size(permutedims(sequence)) == (21, 1)
            @test size(permutedims(annotated_sequence)) == (21, 1)
        end

        @testset "Get residues" begin

            @test getresidues(msa) == M
            @test getresidues(annotated_msa) == M
            @test getresidues(sequence) == S
            @test getresidues(annotated_sequence) == S

            for object in (msa, annotated_msa, sequence, annotated_sequence)
                @test isa(getresidues(object), Matrix{Residue})
            end

            for object in (msa, annotated_msa, sequence, annotated_sequence)
                @test getresiduesequences(msa) == [res"ADEIMSY", res"RCGLFTV", res"NQHKPW-"]
            end
        end

        @testset "Size" begin

            for object in (M, NamedArray(M), msa, annotated_msa)
                @test ncolumns(object) == 7
                @test nsequences(object) == 3
            end

            for object in (S, NamedArray(S), sequence, annotated_sequence)
                @test ncolumns(object) == 21
                @test nsequences(object) == 1
            end
        end

        @testset "Get sequences" begin

            for object in (msa, annotated_msa)
                for seq = 1:3
                    @test getsequence(object, string(seq)) == getsequence(msa, seq)
                    @test size(getsequence(msa, seq)) == (1, 7)
                end
            end
        end

        @testset "Sequence names" begin

            for object in (M, NamedArray(M), msa, annotated_msa)
                @test sequencenames(object) == ["1", "2", "3"]
                # Iterators
                iterator = sequencename_iterator(object)
                @test first(iterator) == "1"
                @test length(iterator) == 3
                @test !isempty(iterator)
                @test collect(iterator) == ["1", "2", "3"]
                @test_throws MethodError iterator[1] # no getindex defined for iterators
            end
        end

        @testset "Rename sequences" begin

            # rename_sequences!
            for object in (NamedArray(M), msa, annotated_msa)
                copied_object = deepcopy(object)
                if isa(copied_object, AnnotatedMultipleSequenceAlignment)
                    setannotsequence!(copied_object, "1", "Name", "One")
                    setannotresidue!(copied_object, "1", "SS", "HHHHHHH")
                end
                new_object = rename_sequences!(copied_object, ["I", "II", "III"])
                @test new_object == copied_object
                @test sequencenames(copied_object) == ["I", "II", "III"]
                if isa(copied_object, AnnotatedMultipleSequenceAlignment)
                    @test getannotsequence(copied_object, "I", "Name") == "One"
                    @test getannotresidue(copied_object, "I", "SS") == "HHHHHHH"
                end
            end

            # rename_sequences
            for object in (NamedArray(M), msa, annotated_msa)
                new_object = rename_sequences(object, ["I", "II", "III"])
                @test sequencenames(object) == ["1", "2", "3"]
                @test sequencenames(new_object) == ["I", "II", "III"]
            end

            # rename one or two sequences
            for object in (NamedArray(M), msa, annotated_msa)
                copied_object = deepcopy(object)
                if isa(copied_object, AnnotatedMultipleSequenceAlignment)
                    setannotsequence!(copied_object, "1", "Name", "One")
                    setannotresidue!(copied_object, "1", "SS", "HHHHHHH")
                end
                new_object = rename_sequences(copied_object, "1" => "I", "2" => "II")
                @test sequencenames(copied_object) == ["1", "2", "3"]
                @test sequencenames(new_object) == ["I", "II", "3"]
                if isa(new_object, AnnotatedMultipleSequenceAlignment)
                    @test getannotsequence(new_object, "I", "Name") == "One"
                    @test getannotresidue(new_object, "I", "SS") == "HHHHHHH"
                end
            end

            # rename_sequences! with Dict
            for object in (NamedArray(M), msa, annotated_msa)
                copied_object = deepcopy(object)
                if isa(copied_object, AnnotatedMultipleSequenceAlignment)
                    setannotsequence!(copied_object, "1", "Name", "One")
                    setannotresidue!(copied_object, "1", "SS", "HHHHHHH")
                end
                rename_sequences!(copied_object, Dict("1" => "I", "2" => "II"))
                @test sequencenames(copied_object) == ["I", "II", "3"]
                if isa(copied_object, AnnotatedMultipleSequenceAlignment)
                    @test getannotsequence(copied_object, "I", "Name") == "One"
                    @test getannotresidue(copied_object, "I", "SS") == "HHHHHHH"
                end
            end

            # rename_sequences! with Pairs
            for object in (NamedArray(M), msa, annotated_msa)
                copied_object = deepcopy(object)
                if isa(copied_object, AnnotatedMultipleSequenceAlignment)
                    setannotsequence!(copied_object, "1", "Name", "One")
                    setannotresidue!(copied_object, "1", "SS", "HHHHHHH")
                end
                rename_sequences!(copied_object, "1" => "I", "2" => "II")
                @test sequencenames(copied_object) == ["I", "II", "3"]
                if isa(copied_object, AnnotatedMultipleSequenceAlignment)
                    @test getannotsequence(copied_object, "I", "Name") == "One"
                    @test getannotresidue(copied_object, "I", "SS") == "HHHHHHH"
                end
            end
        end

        @testset "Rename columns" begin

            # rename_columns!
            for object in (NamedArray(M), msa, annotated_msa)
                copied_object = deepcopy(object)
                oldcols = columnnames(copied_object)
                if isa(copied_object, AnnotatedMultipleSequenceAlignment)
                    setannotfile!(copied_object, "HCat", join(oldcols, ','))
                    hcat_before = getannotfile(copied_object, "HCat")
                end
                new_object =
                    rename_columns!(copied_object, ["A", "B", "C", "D", "E", "F", "G"])
                @test new_object == copied_object
                @test columnnames(copied_object) == ["A", "B", "C", "D", "E", "F", "G"]
                if isa(copied_object, AnnotatedMultipleSequenceAlignment)
                    @test getannotfile(copied_object, "HCat") == hcat_before
                    @test getannotfile(copied_object, "OrigCols") == join(oldcols, ',')
                    @test any(startswith.(keys(getannotfile(copied_object)), "MIToS_"))
                end
            end

            # rename_columns
            for object in (NamedArray(M), msa, annotated_msa)
                if isa(object, AnnotatedMultipleSequenceAlignment)
                    oldcols = columnnames(object)
                    setannotfile!(object, "HCat", join(oldcols, ','))
                end
                new_object = rename_columns(object, ["A", "B", "C", "D", "E", "F", "G"])
                @test columnnames(object) == ["1", "2", "3", "4", "5", "6", "7"]
                @test columnnames(new_object) == ["A", "B", "C", "D", "E", "F", "G"]
                if isa(new_object, AnnotatedMultipleSequenceAlignment)
                    @test getannotfile(new_object, "HCat") ==
                          join(["1", "2", "3", "4", "5", "6", "7"], ',')
                    @test getannotfile(new_object, "OrigCols") ==
                          join(["1", "2", "3", "4", "5", "6", "7"], ',')
                    @test any(startswith.(keys(getannotfile(new_object)), "MIToS_"))
                end
            end

            # rename one or two columns
            for object in (NamedArray(M), msa, annotated_msa)
                copied_object = deepcopy(object)
                if isa(copied_object, AnnotatedMultipleSequenceAlignment)
                    oldcols = columnnames(copied_object)
                    setannotfile!(copied_object, "HCat", join(oldcols, ','))
                    hcat_before = getannotfile(copied_object, "HCat")
                end
                new_object = rename_columns(copied_object, "1" => "A", "2" => "B")
                @test columnnames(copied_object) == ["1", "2", "3", "4", "5", "6", "7"]
                @test columnnames(new_object) == ["A", "B", "3", "4", "5", "6", "7"]
                if isa(new_object, AnnotatedMultipleSequenceAlignment)
                    @test getannotfile(new_object, "HCat") == hcat_before
                    @test getannotfile(new_object, "OrigCols") == join(oldcols, ',')
                    @test any(startswith.(keys(getannotfile(new_object)), "MIToS_"))
                end
            end

            @testset "Renaming after concatenation" begin
                msa_file = joinpath(DATA, "simple.fasta")
                msa_in = read_file(msa_file, FASTA, generatemapping = true)
                h = hcat(msa_in, msa_in)
                before = getannotfile(h, "HCat")
                rename_columns!(h, ["A", "B", "C", "D"])
                @test getannotfile(h, "HCat") == before
                @test getannotfile(h, "OrigCols") == "1_1,1_2,2_1,2_2"
                @test any(startswith.(keys(getannotfile(h)), "MIToS_"))

                v = vcat(h, h)
                before1 = getannotfile(v, "1_HCat")
                before2 = getannotfile(v, "2_HCat")
                rename_columns!(v, ["Q", "W", "E", "R"])
                @test getannotfile(v, "1_HCat") == before1
                @test getannotfile(v, "2_HCat") == before2
                @test getannotfile(v, "OrigCols") == "A,B,C,D"
                @test any(startswith.(keys(getannotfile(v)), "MIToS_"))

                j = join_msas(msa_in, msa_in, [1, 2] .=> [1, 2], axis = 1)
                j_before = getannotfile(j, "HCat")
                rename_columns!(j, ["L", "M", "N", "O"])
                @test getannotfile(j, "HCat") == j_before
                @test getannotfile(j, "OrigCols") == "1_1,1_2,2_1,2_2"
                @test any(startswith.(keys(getannotfile(j)), "MIToS_"))
            end

            # rename_columns! with Dict
            for object in (NamedArray(M), msa, annotated_msa)
                copied_object = deepcopy(object)
                if isa(copied_object, AnnotatedMultipleSequenceAlignment)
                    oldcols = columnnames(copied_object)
                    setannotfile!(copied_object, "HCat", join(oldcols, ','))
                    hcat_before = getannotfile(copied_object, "HCat")
                end
                rename_columns!(copied_object, Dict("1" => "A", "2" => "B"))
                @test columnnames(copied_object) == ["A", "B", "3", "4", "5", "6", "7"]
                if isa(copied_object, AnnotatedMultipleSequenceAlignment)
                    @test getannotfile(copied_object, "HCat") == hcat_before
                    @test getannotfile(copied_object, "OrigCols") == join(oldcols, ',')
                    @test any(startswith.(keys(getannotfile(copied_object)), "MIToS_"))
                end
            end

            # rename_columns! with Pairs
            for object in (NamedArray(M), msa, annotated_msa)
                copied_object = deepcopy(object)
                if isa(copied_object, AnnotatedMultipleSequenceAlignment)
                    oldcols = columnnames(copied_object)
                    setannotfile!(copied_object, "HCat", join(oldcols, ','))
                    hcat_before = getannotfile(copied_object, "HCat")
                end
                rename_columns!(copied_object, "1" => "A", "2" => "B")
                @test columnnames(copied_object) == ["A", "B", "3", "4", "5", "6", "7"]
                if isa(copied_object, AnnotatedMultipleSequenceAlignment)
                    @test getannotfile(copied_object, "HCat") == hcat_before
                    @test getannotfile(copied_object, "OrigCols") == join(oldcols, ',')
                    @test any(startswith.(keys(getannotfile(copied_object)), "MIToS_"))
                end
            end

            @testset "Multiple renames" begin
                obj = deepcopy(annotated_msa)
                oldcols = columnnames(obj)
                before = getannotfile(obj, "HCat")
                rename_columns!(obj, ["A", "B", "C", "D", "E", "F", "G"])
                orig = getannotfile(obj, "OrigCols")
                rename_columns!(obj, ["Q", "W", "E", "R", "T", "Y", "U"])
                @test getannotfile(obj, "OrigCols") == orig
                @test getannotfile(obj, "HCat") == before
            end
        end

        @testset "Column names" begin

            for object in (M, NamedArray(M), msa, annotated_msa)
                @test columnnames(object) == ["1", "2", "3", "4", "5", "6", "7"]
                # Iterators
                iterator = columnname_iterator(object)
                @test first(iterator) == "1"
                @test length(iterator) == 7
                @test !isempty(iterator)
                @test collect(iterator) == ["1", "2", "3", "4", "5", "6", "7"]
                @test_throws MethodError iterator[1] # no getindex defined for iterators
            end
        end

        @testset "Column mapping" begin

            for object in (NamedArray(M), msa, annotated_msa)
                @test getcolumnmapping(object) == [1, 2, 3, 4, 5, 6, 7]
            end
        end

        @testset "Sequence as string" begin

            for object in (M, NamedArray(M), msa, annotated_msa)
                @test stringsequence(object, 1) == "ADEIMSY"
                @test stringsequence(getsequence(object, 1)) == "ADEIMSY"
            end

            for object in (NamedArray(M), msa, annotated_msa)
                @test stringsequence(msa, "1") == "ADEIMSY"
            end
        end

        @testset "Copy and setindex!" begin

            copy_msa = copy(msa)
            deepcopy_msa = deepcopy(msa)

            copy_annotated_msa = copy(annotated_msa)
            deepcopy_annotated_msa = deepcopy(annotated_msa)

            for x in [copy_msa, deepcopy_msa, copy_annotated_msa, deepcopy_annotated_msa]
                x[1, :] = res"YSMIEDA"
                @test vec(x[1, :]) == res"YSMIEDA"

                x["2", :] = res"YSMIEDA"
                @test vec(x["2", :]) == res"YSMIEDA"

                x[:, 1] = res"YYY"
                @test vec(x[:, 1]) == res"YYY"

                x[:, "2"] = res"YYY"
                @test vec(x[:, "2"]) == res"YYY"

                x[end, end] = 'X'
                @test x[end, end] == XAA

                @test length(unique(x)) != 21
            end

            copy_seq = copy(sequence)
            deepcopy_seq = deepcopy(sequence)

            copy_annotated_seq = copy(annotated_sequence)
            deepcopy_annotated_seq = deepcopy(annotated_sequence)

            for x in [copy_seq, deepcopy_seq, copy_annotated_seq, deepcopy_annotated_seq]
                x[1] = XAA
                @test x[1] == XAA

                x["1", "2"] = XAA
                @test x["1", "2"] == XAA

                x[end] = 'X'
                @test x[end] == XAA

                # Special setindex! for sequences:
                x["3"] = GAP
                @test x["3"] == GAP

                @test length(unique(x)) == 19
            end

            @test length(unique(msa)) == 21
            @test length(unique(annotated_msa)) == 21
            @test length(unique(sequence)) == 21
            @test length(unique(annotated_sequence)) == 21
        end
    end

    @testset "MultipleSequenceAlignment Round-trip" begin
        mat = NamedArray(reshape(res"AB", 1, 2))
        msa = MultipleSequenceAlignment(mat)
        @test MultipleSequenceAlignment(msa) === msa
        @test @allocated(MultipleSequenceAlignment(msa)) == 0
    end

    @testset "AnnotatedAlignedSequence Self-identity" begin
        mat = NamedArray(reshape(res"AB", 1, 2))
        ann = Annotations()
        setannotfile!(ann, "Note", "foo")
        aas = AnnotatedAlignedSequence(mat, ann)
        @test AnnotatedAlignedSequence(aas) === aas
        @test annotations(AnnotatedAlignedSequence(aas)) === ann
    end

    @testset "AlignedSequence Self-identity" begin
        mat = NamedArray(reshape(res"AB", 1, 2))
        as = AlignedSequence(mat)
        @test AlignedSequence(as) === as
    end

    @testset "AnnotatedSequence from raw matrix" begin
        rawmat = reshape(res"AC", 1, 2)
        ann = Annotations();
        setannotfile!(ann, "Note", "foo")
        seq_from_raw = AnnotatedSequence(rawmat, ann)
        @test dimnames(seq_from_raw) == ["Seq", "Pos"]
        @test isa(namedmatrix(seq_from_raw), NamedArray)
        @test annotations(seq_from_raw) === ann # it wraps the annotation object
        @test getresidues(seq_from_raw) == rawmat
    end

    @testset "AnnotatedSequence from AlignedSequence" begin
        mat_gap = NamedArray(reshape(res"A-C", 1, 3))
        alnseq = AlignedSequence(mat_gap)
        seq_from_aln = AnnotatedSequence(alnseq)
        @test dimnames(seq_from_aln) == ["Seq", "Pos"]
        @test isempty(annotations(seq_from_aln))
        @test all(getresidues(seq_from_aln) .!= GAP)
    end

    @testset "AnnotatedSequence from AnnotatedAlignedSequence" begin
        mat_gap = NamedArray(reshape(res"A-C", 1, 3))
        ann = Annotations();
        setannotfile!(ann, "Note", "foo")
        aas = AnnotatedAlignedSequence(mat_gap, ann) # has gaps as it is aligned
        seq_from_aas = AnnotatedSequence(aas) # deletes gaps to get the unaligned sequence
        @test annotations(seq_from_aas) == ann && annotations(seq_from_aas) !== ann
    end

    @testset "AnnotatedSequence no-op identity" begin
        seq = AnnotatedSequence(res"AC")
        orig_mat = namedmatrix(seq)
        orig_ann = annotations(seq)
        @test AnnotatedSequence(seq) === seq
        @test namedmatrix(seq) === orig_mat && annotations(seq) === orig_ann
    end

    @testset "Branch coverage inside AnnotatedSequence(matrix,annot)" begin
        named_col = NamedArray(reshape(res"AD", 1, 2), (["id"], ["1", "2"]), ("Seq", "Col"))
        seq_branch1 = AnnotatedSequence(named_col, Annotations())
        @test dimnames(named_col)[2] == "Col"
        @test dimnames(seq_branch1)[2] == "Pos"


        named_pos = NamedArray(reshape(res"AD", 1, 2), (["id"], ["1", "2"]), ("Seq", "Pos"))
        seq_branch2 = AnnotatedSequence(named_pos, Annotations())
        @test dimnames(named_pos)[2] == "Pos"
        @test dimnames(seq_branch2)[2] == "Pos"
    end
end
