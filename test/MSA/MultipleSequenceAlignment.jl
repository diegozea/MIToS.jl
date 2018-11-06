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

        M = reshape(reinterpret(Residue,collect(1:21)),(3,7))
        m = NamedArray(M)

        T = permutedims(M, [2,1])

        S = reshape(reinterpret(Residue,collect(1:21)),(1,21))
        s = NamedArray(S)

        @test MultipleSequenceAlignment(M) == MultipleSequenceAlignment(m)
        @test AnnotatedMultipleSequenceAlignment(M) == AnnotatedMultipleSequenceAlignment(m)

        @test AlignedSequence(S) == AlignedSequence(s)
        @test AnnotatedAlignedSequence(S) == AnnotatedAlignedSequence(s)

        @test MultipleSequenceAlignment(M) != MultipleSequenceAlignment(T)
        @test AnnotatedMultipleSequenceAlignment(M) != AnnotatedMultipleSequenceAlignment(T)

        @test_throws AssertionError AlignedSequence(M)
        @test_throws AssertionError AnnotatedAlignedSequence(M)
    end

    @testset "MSA & sequences" begin

        M = reshape(reinterpret(Residue,collect(1:21)),(3,7))
        msa = MultipleSequenceAlignment(M)
        annotated_msa = AnnotatedMultipleSequenceAlignment(M)

        S = reshape(reinterpret(Residue,collect(1:21)),(1,21))
        sequence = AlignedSequence(S)
        annotated_sequence = AnnotatedAlignedSequence(S)

        @testset "Getters" begin

            @test isa(annotations(annotated_msa), Annotations)
            @test isa(annotations(annotated_sequence), Annotations)

            for object in (msa, annotated_msa, sequence, annotated_sequence)
                @test isa(namedmatrix(object), NamedResidueMatrix{Array{Residue,2}})
            end
        end

        @testset "Dimension names" begin

            for object in (msa, annotated_msa, sequence, annotated_sequence)
                @test dimnames(namedmatrix(object)) == ["Seq","Col"]
            end
        end

        @testset "Convert" begin

            msa2annot = convert(AnnotatedMultipleSequenceAlignment, msa)
            annot2msa = convert(MultipleSequenceAlignment, annotated_msa)

            seq2annot = convert(AnnotatedAlignedSequence, sequence)
            annot2seq = convert(AlignedSequence, annotated_sequence)

            @test isa(annotations(msa2annot), Annotations)
            @test isa(annotations(seq2annot), Annotations)

            @test_throws MethodError annotations(annot2msa)
            @test_throws MethodError annotations(annot2seq)

            @test namedmatrix(msa2annot) == namedmatrix(msa)
            @test namedmatrix(annot2msa) == namedmatrix(annotated_msa)
            @test namedmatrix(seq2annot) == namedmatrix(sequence)
            @test namedmatrix(annot2seq) == namedmatrix(annotated_sequence)
        end

        @testset "AbstractArray Interface" begin

            @test size(msa) == (3,7)
            @test size(annotated_msa) == (3,7)

            @test size(sequence) == (1,21)
            @test size(annotated_sequence) == (1,21)

            for object in (msa, annotated_msa, sequence, annotated_sequence)
                @test length(object) == 21
            end
        end

        @testset "Indexing" begin

            for object in (msa, annotated_msa)
                for i in 1:3, j in 1:7
                    @test object[string(i),string(j)] == object[i,j]
                end
            end

            for object in (sequence, annotated_sequence)
                for j in 1:7
                    @test object["1",string(j)] == object[1,j]
                    @test object[j] == object[1,j]
                    # Special sequence indexing:
                    @test object[string(j)] == object[1,j]
                end
            end
        end

        @testset "Show" begin

            out = IOBuffer()

            show(out, msa)
            str = String(take!(out))
            @test startswith(str,"MultipleSequenceAlignment : ")
            @test occursin("Seq", str)
            @test occursin("Col", str)
            @test length(split(str,'\n')) == 6

            show(out, annotated_msa)
            str = String(take!(out))
            @test startswith(str,"AnnotatedMultipleSequenceAlignment with 0 annotations : ")
            @test occursin("Seq", str)
            @test occursin("Col", str)
            @test length(split(str,'\n')) == 6

            show(out, sequence)
            str = String(take!(out))
            @test startswith(str,"AlignedSequence : ")
            @test occursin("Seq", str)
            @test occursin("Col", str)
            @test length(split(str,'\n')) == 4

            show(out, annotated_sequence)
            str = String(take!(out))
            @test startswith(str,"AnnotatedAlignedSequence with 0 annotations : ")
            @test occursin("Seq", str)
            @test occursin("Col", str)
            @test length(split(str,'\n')) == 4
        end

        @testset "Transpose" begin

            @test transpose(msa) == permutedims(msa)
            @test size(transpose(msa)) == (7,3)
            @test size(transpose(annotated_msa)) == (7,3)

            @test transpose(sequence) == permutedims(sequence)
            @test size(transpose(sequence)) == (21,1)
            @test size(transpose(annotated_sequence)) == (21,1)
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
                @test getresiduesequences(msa) == [res"ADEIMSY",
                                                   res"RCGLFTV",
                                                   res"NQHKPW-"]
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
                for seq in 1:3
                    @test getsequence(object, string(seq)) == getsequence(msa, seq)
                    @test size(getsequence(msa, seq)) == (1,7)
                end
            end
        end

        @testset "Sequence names" begin

            for object in (M, NamedArray(M), msa, annotated_msa)
                @test sequencenames(object) == ["1","2","3"]
            end
        end

        @testset "Column names" begin

            for object in (M, NamedArray(M), msa, annotated_msa)
                @test columnnames(object) == ["1","2","3","4","5","6","7"]
            end
        end

        @testset "Column mapping" begin

            for object in (NamedArray(M), msa, annotated_msa)
                @test getcolumnmapping(object) == [1,2,3,4,5,6,7]
            end
        end

        @testset "Sequence as string" begin

            for object in (M, NamedArray(M), msa, annotated_msa)
                @test stringsequence(object, 1) == "ADEIMSY"
                @test stringsequence(getsequence(object,1)) == "ADEIMSY"
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
                x[1,:] = res"YSMIEDA"
                @test vec(x[1,:]) == res"YSMIEDA"

                x["2",:] = res"YSMIEDA"
                @test vec(x["2",:]) == res"YSMIEDA"

                x[:,1] = res"YYY"
                @test vec(x[:,1]) == res"YYY"

                x[:,"2"] = res"YYY"
                @test vec(x[:,"2"]) == res"YYY"

                x[end,end] = 'X'
                @test x[end,end] == XAA

                @test length(unique(x)) != 21
            end

            copy_seq = copy(sequence)
            deepcopy_seq = deepcopy(sequence)

            copy_annotated_seq = copy(annotated_sequence)
            deepcopy_annotated_seq = deepcopy(annotated_sequence)

            for x in [copy_seq, deepcopy_seq, copy_annotated_seq, deepcopy_annotated_seq]
                x[1] = XAA
                @test x[1] == XAA

                x["1","2"] = XAA
                @test x["1","2"] == XAA

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
end
