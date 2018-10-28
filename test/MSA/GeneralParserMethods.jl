@testset "GeneralParserMethods" begin

    @testset "To parse MSA & mapping" begin

        sequences = ["ADEIMSY","RCGLFTV","NQHKPW-"]
        matrixres = reshape(reinterpret(Residue, collect(1:21)), (3,7))
        msa, map = MSA._to_msa_mapping(sequences)

        @test getarray(msa) == matrixres
        @test map == ["1,2,3,4,5,6,7", "1,2,3,4,5,6,7", "1,2,3,4,5,6,"]
        @test sequencenames(msa) == ["1", "2", "3"]
        # MSA constructors adds dimension names
        # @test dimnames(msa) == [:Seq, :Col]

        msa, map = MSA._to_msa_mapping(sequences,["a/11-17","b/11-17","c/11-16"])

        @test getarray(msa) == matrixres
        @test map == ["11,12,13,14,15,16,17", "11,12,13,14,15,16,17", "11,12,13,14,15,16,"]
        @test sequencenames(msa) == ["a/11-17", "b/11-17", "c/11-16"]
        # MSA constructors adds dimension names
        # @test dimnames(msa) == [:Seq, :Col]

        @test_throws ErrorException MSA._to_msa_mapping(sequences,["a/1-7","b/1-7","c/1-7"])
        @test_throws ErrorException MSA._to_msa_mapping(["AD","EI","MSY"])
    end

    @testset "Delete full gap columns" begin

        M = reshape(reinterpret(Residue,collect(1:21)),(3,7))
        M[:,[2,4,6]] .= GAP
        msa = MultipleSequenceAlignment(M)
        named = namedmatrix(msa)
        annotated_msa = AnnotatedMultipleSequenceAlignment(M)

        d_M = deletefullgapcolumns(M)
        d_msa = deletefullgapcolumns(msa)
        d_named = deletefullgapcolumns(named)
        d_annotated_msa = deletefullgapcolumns(annotated_msa)

        @test size(d_M) == (3,4)

        for object in [d_msa, d_named, d_annotated_msa]
            @test size(object) == (3,4)
            @test getcolumnmapping(object) == [1,3,5,7]
        end

        d_msa = deletefullgapcolumns!(msa)
        d_annotated_msa = deletefullgapcolumns!(annotated_msa)

        @test d_msa == msa
        @test d_annotated_msa == annotated_msa
    end
end
