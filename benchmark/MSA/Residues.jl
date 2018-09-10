# En Pfam 30.0 PF00400 has 268378 sequences
@benchgroup "Residue conversions" ["IO", "MSA"] begin

    Random.seed!(1)
    chars = rand(['.','-','a':'z'...,'A':'Z'...], 268378*2)
    residues = Residue[ char for char in chars ]
    ints = Int[ res for res in residues ]

    @bench "char2res" begin
                        for char in $chars
                            Residue(char)
                        end
                      end
    @bench "res2char" begin
                        for res in $residues
                            Char(res)
                        end
                      end
    @bench "int2res"  begin
                        for i in $ints
                            Residue(i)
                        end
                      end
    @bench "res2int"  begin
                        for res in $residues
                            Int(res)
                        end
                      end
end
