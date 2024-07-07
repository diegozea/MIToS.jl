let ascii_str = "#=GR O31698/18-71 SS    CCCHHHHHHHHHHHHHHHEEEEEEEEEEEEEEEEHHH\n",
    utf8_str = "#=GF CC   (Römling U.  and Galperin M.Y. “Bacterial cellulose\n"

    SUITE["Utils"]["get_n_words"]["ascii"] = @benchmarkable get_n_words($ascii_str, 3)
    SUITE["Utils"]["get_n_words"]["utf8"] = @benchmarkable get_n_words($utf8_str, 3)
end
