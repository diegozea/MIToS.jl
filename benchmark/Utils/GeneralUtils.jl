let ascii_str = "#=GR O31698/18-71 SS    CCCHHHHHHHHHHHHHHHEEEEEEEEEEEEEEEEHHH\n",
    utf8_str = "#=GF CC   (Römling U.  and Galperin M.Y. “Bacterial cellulose\n"

    SUITE["Utils"]["get_n_words"]["ascii"] = @benchmarkable get_n_words($ascii_str, 3)
    SUITE["Utils"]["get_n_words"]["utf8"] = @benchmarkable get_n_words($utf8_str, 3)
end

begin
    SUITE["Utils"]["hascoordinates"]["valid"] =
        @benchmarkable hascoordinates("O83071/192-246")
    SUITE["Utils"]["hascoordinates"]["invalid"] =
        @benchmarkable hascoordinates("invalid_format")
end

let mat = rand(500, 500), side = 500, vec = matrix2list(mat)
    SUITE["Utils"]["matrix2list"]["upper"] = @benchmarkable matrix2list($mat)
    SUITE["Utils"]["matrix2list"]["upper_diagonal"] =
        @benchmarkable matrix2list($mat; diagonal = true)
    SUITE["Utils"]["list2matrix"]["upper"] = @benchmarkable list2matrix($vec, $side)
    SUITE["Utils"]["list2matrix"]["upper_diagonal"] = @benchmarkable list2matrix(
        matrix2list($mat; diagonal = true),
        $side;
        diagonal = true,
    )
end