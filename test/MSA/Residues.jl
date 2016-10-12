@testset "Residue" begin

    @testset "Comparisons" begin

        @test Residue(1) == Residue(1)

        @test GAP == GAP
        @test XAA == XAA
        @test INV == INV

        @test GAP != INV
        @test GAP != XAA

        for i in 1:23, j in 1:23
            if j != i
                @test Residue(i) != Residue(j)
            else
                @test Residue(i) == Residue(j)
            end
        end
    end

    @testset "Convert" begin

        alphabet = "ARNDCQEGHILKMFPSTWYV"

        @test Int[ Residue(char) for char in alphabet] == Int[ i for i in 1:20 ]
        @test Int[ Residue(char) for char in lowercase(alphabet)] == Int[ 21 for i in 1:20 ]

        @testset "GAP" begin

            @test Int(GAP) == 21

            @test Char(GAP) == '-'

            @test GAP == Residue('-')
            @test GAP == Residue('*')
            @test GAP == Residue('.')

            for X in "UOBZJX"
                @test GAP != Residue(X)
                @test GAP == Residue(lowercase(X))
            end

            @test GAP != Residue(42)
            @test GAP != Residue(-1)
        end

        @testset "XAA" begin

            @test Int(XAA) == 22

            @test Char(XAA) == 'X'

            for X in "UOBZJX"
                @test XAA == Residue(X)
                @test XAA != Residue(lowercase(X))
            end
        end

        @testset "INV" begin

            @test Int(INV) == 23

            @test Char(INV) == '\Ufffd'

            @test INV == Residue('\Ufffd')
            @test INV == Residue('ñ')
            @test INV == Residue(42)
            @test INV == Residue(-1)

            for X in "UOBZJX"
                @test INV != Residue(X)
                @test INV != Residue(lowercase(X))
            end
        end
    end

    @testset "Valid Residues" begin

        for res in res"ARNDCQEGHILKMFPSTWYV-X"
            @test isvalid(res)
        end

        @test !isvalid(INV)
        @test !isvalid(Residue('ñ'))
        @test !isvalid(reinterpret(Residue, 42))
    end

    @testset "Strings" begin

        alphabet = "ARNDCQEGHILKMFPSTWYV-X�"

        @test res"ARNDCQEGHILKMFPSTWYV-X�" == Residue[ char for char in alphabet ]
        @test String(res"ARNDCQEGHILKMFPSTWYV-X�") == alphabet
    end

    @testset "String vector" begin

        msa =  ["DAWAEF",
                "DAWAED",
                "DAYCMD" ]

        badmsa=["DAWAEF",
                "DAWAED",
                "DAM"    ]

        @test convert(Matrix{Residue}, msa) == Residue['D' 'A' 'W' 'A' 'E' 'F'
                                                       'D' 'A' 'W' 'A' 'E' 'D'
                                                       'D' 'A' 'Y' 'C' 'M' 'D']

        @test_throws ErrorException convert(Matrix{Residue}, String[])
        @test_throws ErrorException convert(Matrix{Residue}, badmsa)
    end

    @testset "Random" begin

        for i in 1:5000
            res = rand(Residue)
            @test res != GAP
            @test res != XAA
            @test res != INV
            @test isvalid(res)
            @test typeof(res) == Residue
        end
    end

    @testset "Initialized arrays" begin

        @test size(rand(Residue, 20)) == (20,)
        @test typeof(rand(Residue, 20)) == Array{Residue, 1}
        @test size(rand(Residue, 20,30)) == (20,30)
        @test typeof(rand(Residue, 20, 30)) == Array{Residue, 2}

        @test zero(Residue) == GAP
        @test one( Residue) == XAA

        @test zeros(Residue,4) == [GAP, GAP, GAP, GAP]
        @test ones(Residue,4)  == [XAA, XAA, XAA, XAA]

        @test zeros(Residue,2,2) == [   GAP GAP
                                        GAP GAP ]

        @test ones(Residue,2,2)  == [   XAA XAA
                                        XAA XAA ]
    end

    @testset "Other base methods" begin

        for i in 1:23
            @test bits(i) == bits(Residue(i))
        end

        @test typemin(Residue) == Residue(1)
        @test typemax(Residue) == INV
    end
end
