@benchgroup "get_n_words" ["IO", "Pfam"] begin

    ascii_str = "#=GR O31698/18-71 SS    CCCHHHHHHHHHHHHHHHEEEEEEEEEEEEEEEEHHH\n"
    utf8_str  = "#=GF CC   (Römling U.  and Galperin M.Y. “Bacterial cellulose\n"

    @bench "ascii" get_n_words($ascii_str, 3)
    @bench "utf8"  get_n_words($utf8_str,  3)
end
