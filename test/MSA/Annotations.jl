@testset "Annotations" begin

    @testset "Empty annotations" begin

        annot = Annotations()

        @test length(annot.file) == 0
        @test length(annot.sequences) == 0
        @test length(annot.columns) == 0
        @test length(annot.residues) == 0

        @test length(annot) == 0

        @test ncolumns(annot) == -1

        @test isempty(annot)
    end

    @testset "Getters & Setters" begin

        annot = Annotations()
        example_str = "CCCHHHHHHHHHHHHHHHEEEEEEEEEEEEEEEEHHH"

        setannotresidue!(annot, "O31698/18-71", "SS", example_str)

        @test_throws AssertionError setannotresidue!(annot, "O31698/18-71",
                                                    String(rand('A':'Z',51)), example_str)

        @test_throws AssertionError setannotresidue!(annot, "O31698/18-71", "Feature Name",
                                                    example_str)

        @test ncolumns(annot) == 37

        setannotfile!(annot, "AC", "PF00571")
        setannotcolumn!(annot, "SS_cons", example_str)
        setannotsequence!(annot, "O31698/88-139", "OS", "Bacillus subtilis")

        @test getannotfile(annot, "AC") == "PF00571"
        @test getannotcolumn(annot, "SS_cons") == example_str
        @test getannotsequence(annot, "O31698/88-139", "OS") == "Bacillus subtilis"
        @test getannotresidue(annot, "O31698/18-71", "SS") == example_str

        @test getannotfile(annot, "An", "Default") == "Default"
        @test getannotcolumn(annot, "Other", "Default") == "Default"
        @test getannotsequence(annot, "O31698/1-88", "OS", "Default") == "Default"
        @test getannotresidue(annot, "O31698/1-88", "SS", "Default") == "Default"

        @test ncolumns(annot) == 37

        @test_throws DimensionMismatch  setannotresidue!(annot,
                                        "O31698/18-71", "AS", "__*__")

        @test_throws DimensionMismatch  setannotcolumn!(annot, "SS_cons",
                                        "---CCCCCHHHHHHHHHHHHHEEEEEEEEEEEEEEEEEEH---")
    end

    @testset "Copy, deepcopy and empty!" begin

        annot = Annotations()
        setannotresidue!(annot,"O31698/18-71","SS","CCCHHHHHHHHHHHHHHHEEEEEEEEEEEEEEEEHHH")
        setannotfile!(annot, "AC", "PF00571")
        setannotcolumn!(annot, "SS_cons", "CCCCCHHHHHHHHHHHHHEEEEEEEEEEEEEEEEEEH")
        setannotsequence!(annot, "O31698/88-139", "OS", "Bacillus subtilis")

        copy_annot = copy(annot)
        @test copy_annot == annot
        empty!(copy_annot)
        @test ncolumns(annot) == 37
        @test ncolumns(copy_annot) == -1

        @test !isempty(annot)
        @test isempty(copy_annot)

        deepcopy_annot = deepcopy(annot)
        @test deepcopy_annot == annot
        empty!(deepcopy_annot)
        @test ncolumns(annot) == 37
        @test ncolumns(deepcopy_annot) == -1

        @test !isempty(annot)
        @test isempty(deepcopy_annot)
    end

    @testset "Filter" begin

        annot = Annotations()
        setannotresidue!(annot,"O31698/18-71","SS","CCCHHHHHHHHHHHHHHHEEEEEEEEEEEEEEEEHHH")
        setannotfile!(annot, "AC", "PF00571")
        setannotcolumn!(annot, "SS_cons", "CCCCCHHHHHHHHHHHHHEEEEEEEEEEEEEEEEEEH")
        setannotsequence!(annot, "O31698/88-139", "OS", "Bacillus subtilis")

        @test_throws AssertionError filtercolumns!(annot, [true, false, true])

        #filtersequences!(annot, IndexedArray(["O31698/88-139", "O31698/18-71"]), [false, true])
        #@test length( getannotsequence(annot) ) == 0
        #filtersequences!(annot, IndexedArray(["O31698/88-139", "O31698/18-71"]), [true, false])
        #@test length( getannotresidue(annot) ) == 0

        mask = collect("CCCCCHHHHHHHHHHHHHEEEEEEEEEEEEEEEEEEH") .!= Ref('E')
        filtercolumns!(annot, mask)
        @test ncolumns(annot) == 19
        @test getannotcolumn(annot, "SS_cons") == "CCCCCHHHHHHHHHHHHHH"

        filtercolumns!(annot, [1,2,19])
        @test ncolumns(annot) == 3
        @test getannotcolumn(annot, "SS_cons") == "CCH"
    end

    @testset "Print" begin

        output_string = """
                        #=GF AC	PF00571
                        #=GS O31698/88-139	OS	Bacillus subtilis
                        #=GR O31698/18-71	SS	CCCHHHHHHHHHHHHHHHEEEEEEEEEEEEEEEEHHH
                        #=GC SS_cons			CCCCCHHHHHHHHHHHHHEEEEEEEEEEEEEEEEEEH
                        """

        io = IOBuffer()

        annot = Annotations()
        setannotresidue!(annot,"O31698/18-71","SS","CCCHHHHHHHHHHHHHHHEEEEEEEEEEEEEEEEHHH")
        setannotfile!(annot, "AC", "PF00571")
        setannotcolumn!(annot, "SS_cons", "CCCCCHHHHHHHHHHHHHEEEEEEEEEEEEEEEEEEH")
        setannotsequence!(annot, "O31698/88-139", "OS", "Bacillus subtilis")

        print(io, annot)
        @test String(take!(io)) == output_string

        show(io, annot)
        @test String(take!(io)) == output_string
    end
end
