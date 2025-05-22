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

@testset "_get_function_name" begin
    @test MIToS.Utils._get_function_name("Module.SubModule.function_name") == "function_name"
    @test MIToS.Utils._get_function_name("function_name") == "function_name"
    @test MIToS.Utils._get_function_name("Module.function_name") == "function_name"
end

@testset "hascoordinates" begin

    @test hascoordinates("O83071/192-246")
    @test !hascoordinates("O83071")
end

@testset "select_element" begin

    @test select_element([1]) == 1
    @test_throws ErrorException select_element([])
    @test_logs (:warn, r"There are more than one") select_element([1,2])
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
        @test !isnotemptyfile("nonexistentfile") # test for isnotemptyfile with non existent file

        @test check_file(file_path) == file_path
        @test_throws ErrorException check_file("nonexistentfile")
        # Test @warn message for empty file
        empty_file_path = joinpath(DATA, "empty_test_file.txt")
        touch(empty_file_path)
        try
            @test_logs (:warn, r"is empty!") check_file(empty_file_path)
        finally
            rm(empty_file_path)
        end
    end

    # Mock FileFormat for testing read_file and write_file
    struct DummyFormat <: MIToS.Utils.FileFormat end
    function MIToS.Utils.parse_file(io::IO, ::Type{DummyFormat})
        return read(io, String)
    end
    function MIToS.Utils.print_file(io::IO, data::String, ::Type{DummyFormat})
        print(io, data)
    end

    @testset "read_file" begin
        # Test reading from a local file
        dummy_content = "This is a dummy file."
        dummy_file_path = "dummy_test_file.txt"
        open(dummy_file_path, "w") do f
            write(f, dummy_content)
        end
        try
            @test read_file(dummy_file_path, DummyFormat) == dummy_content
        finally
            rm(dummy_file_path)
        end

        # Test reading a gzipped local file
        gzipped_dummy_file_path = "dummy_test_file.txt.gz"
        gz_stream = GzipCompressorStream(open(gzipped_dummy_file_path, "w"))
        write(gz_stream, dummy_content)
        close(gz_stream)
        try
            @test read_file(gzipped_dummy_file_path, DummyFormat) == dummy_content
        finally
            rm(gzipped_dummy_file_path)
        end

        # Test reading from a URL (requires internet)
        # This test uses a reliable small public file.
        # If this is problematic, this specific test might need to be skipped or mocked.
        # For now, assuming internet access based on existing download tests.
        # Using a UniProt FASTA file as an example, expecting it to be parsed as String by DummyFormat
        # This will not actually parse FASTA, just read the content as a string.
        # The test is more about the mechanism of reading from URL.
        # A small, stable text file URL would be better if available.
        # Using httpbin to get a small text response
        httpbin_url = "http://httpbin.org/robots.txt"
        try
            # The content of robots.txt is known and stable
            expected_robots_txt_content = "User-agent: *\nDisallow: /deny\n"
            @test read_file(httpbin_url, DummyFormat) == expected_robots_txt_content
        catch e
            @warn "URL read_file test failed. It might be due to network issues or the resource being unavailable." e
        end
    end

    @testset "write_file" begin
        dummy_content = "This is content for write_file."
        dummy_file_path = "write_test_file.txt"
        gzipped_dummy_file_path = "write_test_file.txt.gz"

        # Test writing to a local file
        try
            write_file(dummy_file_path, dummy_content, DummyFormat)
            @test read(dummy_file_path, String) == dummy_content
        finally
            isfile(dummy_file_path) && rm(dummy_file_path)
        end

        # Test writing to a gzipped local file
        try
            write_file(gzipped_dummy_file_path, dummy_content, DummyFormat)
            # Verify by reading back through GzipDecompressorStream
            read_content = open(gzipped_dummy_file_path, "r") do f
                gz_stream = GzipDecompressorStream(f)
                read(gz_stream, String)
            end
            @test read_content == dummy_content
            # Also check if it's a valid gzip file using _check_gzip_file
            @test MIToS.Utils._check_gzip_file(gzipped_dummy_file_path) == gzipped_dummy_file_path
        finally
            isfile(gzipped_dummy_file_path) && rm(gzipped_dummy_file_path)
        end
    end

    @testset "Deprecation warnings for write/print" begin
        dummy_content = "Deprecated write/print test."
        dummy_file_path = "deprecated_write_test.txt"

        # Test deprecated write
        try
            @test_deprecated MIToS.Utils.write(dummy_file_path, dummy_content, DummyFormat)
            @test read(dummy_file_path, String) == dummy_content
        finally
            isfile(dummy_file_path) && rm(dummy_file_path)
        end

        # Test deprecated print
        io = IOBuffer()
        @test_deprecated MIToS.Utils.print(io, dummy_content, DummyFormat)
        @test String(take!(io)) == dummy_content
        # Test deprecated print to stdout (implicitly)
        # This is harder to capture without interfering with test runner output
        # We rely on the previous IOBuffer test to cover the print deprecation logic
    end

    @testset "Deprecation warning for deleteitems!" begin
        # From src/Utils/Utils.jl
        # function deleteitems!(vector::Vector, items)
        #     Base.depwarn("`deleteitems!` is deprecated, use `filter!(x -> !(x in items), vector)` instead.", :deleteitems!)
        #     filter!(x -> !(x in items), vector)
        # end
        # Test is directly here as it's a small, self-contained deprecated function from Utils.jl
        # No separate Utils.jl test file exists, and this is part of the GeneralUtils tests context
        v = [1, 2, 3, 4, 5]
        items_to_delete = [2, 4]
        @test_deprecated MIToS.Utils.deleteitems!(v, items_to_delete)
        @test v == [1, 3, 5]

        v_str = ["a", "b", "c", "d"]
        items_to_delete_str = ["b", "d"]
        @test_deprecated MIToS.Utils.deleteitems!(v_str, items_to_delete_str)
        @test v_str == ["a", "c"]
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




