# En Pfam 30.0 PF00400 has 268378 sequences
@benchgroup "Residue conversions" ["IO", "MSA"] begin

    srand(1)
    chars = rand(['.','-','a':'z'...,'A':'Z'...], 268378)
    residues = Residue[ char for char in chars ]
    ints = Int[ res for res in residues ]

    @bench "char2res" foreach(Residue, $chars)
    @bench "res2char" foreach(Char, $residues)
    @bench "int2res"  foreach(Residue, $ints)
    @bench "res2int"  foreach(Int,$residues)
end
