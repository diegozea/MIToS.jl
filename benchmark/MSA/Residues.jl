# En Pfam 30.0 PF00400 has 268378 sequences
let chars = rand(Random.MersenneTwister(1), ['.', '-', 'a':'z'..., 'A':'Z'...], 268378 * 2),
    residues = Residue[char for char in chars],
    ints = Int[res for res in residues]

    SUITE["MSA"]["Residue conversions"]["char2res"] = @benchmarkable Residue.($chars)
    SUITE["MSA"]["Residue conversions"]["res2char"] = @benchmarkable Char.($residues)
    SUITE["MSA"]["Residue conversions"]["int2res"] = @benchmarkable Residue.($ints)
    SUITE["MSA"]["Residue conversions"]["res2int"] = @benchmarkable Int.($residues)
end
