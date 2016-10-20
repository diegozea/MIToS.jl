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

    @testset "Creation and conversion" begin

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
    end


end
