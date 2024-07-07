@testset "get_n_words!" begin

    line = "#=GF AC PF00571"
    @test get_n_words(line, 1) == String[line]
    @test get_n_words(line, 2) == String["#=GF", "AC PF00571"]
    @test get_n_words(line, 3) == String["#=GF", "AC", "PF00571"]
    @test get_n_words(line, 4) == String["#=GF", "AC", "PF00571"]

    @test get_n_words("\n", 1) == String["\n"]
    @test get_n_words("#", 1) == String["#"]

    # ASCII
    str = "#=GR O31698/18-71 SS    CCCHHHHHHHHHHHHHHHEEEEEEEEEEEEEEEEHHH"
    @test get_n_words(str, 3) ==
          String["#=GR", "O31698/18-71", "SS    CCCHHHHHHHHHHHHHHHEEEEEEEEEEEEEEEEHHH"]

    # UTF-8
    str = "#=GF CC   (Römling U.  and Galperin M.Y. “Bacterial cellulose"
    @test get_n_words(str, 3) ==
          String["#=GF", "CC", "(Römling U.  and Galperin M.Y. “Bacterial cellulose"]

    str = "#=GF CC   not present in all SecA2–SecY2 systems. This family of Asp5 is"
    @test get_n_words(str, 3) == String[
        "#=GF",
        "CC",
        "not present in all SecA2–SecY2 systems. This family of Asp5 is",
    ]
end

@testset "hascoordinates" begin

    @test hascoordinates("O83071/192-246")
    @test !hascoordinates("O83071")
end

@testset "select_element" begin

    @test select_element([1]) == 1
    @test_throws ErrorException select_element([])
end

@testset "Matrices to and from lists" begin

    @testset "matrix2list" begin

        mat = [
            1 2 3
            4 5 6
            7 8 9
        ]

        @test matrix2list(mat) == [2, 3, 6]
        @test matrix2list(mat, diagonal = true) == [1, 2, 3, 5, 6, 9]
        @test matrix2list(mat, part = "lower") == [4, 7, 8]
        @test matrix2list(mat, part = "lower", diagonal = true) == [1, 4, 7, 5, 8, 9]
    end

    @testset "list2matrix" begin

        mat = [
            1 2 3
            2 5 6
            3 6 9
        ]

        @test triu(list2matrix([2, 3, 6], 3), 1) == triu(mat, 1)
        @test list2matrix([1, 2, 3, 5, 6, 9], 3, diagonal = true) == mat
    end
end

@testset "General IO" begin

    @testset "lineiterator" begin

        # Julia 0.6: eachline return lines without line endings by default
        ppap = "pen\npineapple\napple\npen\n"
        @test collect(lineiterator(ppap)) == collect(eachline(IOBuffer(ppap)))

        @test collect(lineiterator("Hola")) == ["Hola"]
        @test collect(lineiterator("Hola\n")) == ["Hola"]
        @test collect(lineiterator("\n")) == [""]
        @test collect(lineiterator("Hola\nMundo")) == ["Hola", "Mundo"]
        @test collect(lineiterator("Hola\nMundo\n")) == ["Hola", "Mundo"]
        @test collect(lineiterator("Hola\nMundo\n\n")) == ["Hola", "Mundo", ""]
    end

    @testset "File checking" begin

        file_path = joinpath(DATA, "simple.fasta")

        @test isnotemptyfile(file_path)
        @test !isnotemptyfile(joinpath(DATA, "emptyfile"))

        @test check_file(file_path) == file_path
        @test_throws ErrorException check_file("nonexistentfile")
    end

    @testset "Download file" begin

        try
            @test ".tmp" == download_file(
                "http://www.uniprot.org/uniprot/P69905.fasta",
                ".tmp",
                headers = Dict(
                    "User-Agent" => "Mozilla/5.0 (compatible; MSIE 7.01; Windows NT 5.0)",
                ),
            )
        finally
            if isfile(".tmp")
                rm(".tmp")
            end
        end
    end

    @testset "Download file: redirect" begin
        try
            # Use https://httpbin.io/ and example.com to test redirection
            download_file("https://httpbin.io/redirect-to?url=https://example.com", ".tmp")
            @test occursin("Example Domain", read(".tmp", String))
        finally
            isfile(".tmp") && rm(".tmp")
        end
    end

    @testset "Test _check_gzip_file" begin
        for file in readdir(DATA)
            filename = joinpath(DATA, file)
            if file != "2vqc.xml.gz" # is a decompressed file that has a wrong .gz extension
                @test MIToS.Utils._check_gzip_file(filename) == filename
            else
                @test_throws ErrorException MIToS.Utils._check_gzip_file(filename)
            end
        end
    end

    @testset "Download a gz file" begin
        # Use https://www.rcsb.org/pdb/files/3NIR.pdb.gz to test downloading a gz file
        # without a filename
        filename = ""
        try
            filename = download_file("https://www.rcsb.org/pdb/files/3NIR.pdb.gz")
            @test endswith(filename, ".gz")
            @test MIToS.Utils._check_gzip_file(filename) == filename
        finally
            if isfile(filename)
                rm(filename)
            end
        end
    end
end




