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
                @test isa(namedmatrix(object), NamedArray{Residue,2})
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
    end

end
