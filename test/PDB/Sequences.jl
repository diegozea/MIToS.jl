@testset "Check amino acid residues" begin
    @test PDB._is_aminoacid("MET")  # methionine
    @test PDB._is_aminoacid("MSE")  # selenomethionine
    @test !PDB._is_aminoacid("HOH") # water
    @test !PDB._is_aminoacid("SO4") # sulfate
    @test !PDB._is_aminoacid("MG")  # magnesium
    @test !PDB._is_aminoacid("CA")  # calcium
    @test !PDB._is_aminoacid("A")   # adenine
    @test PDB._is_aminoacid("ALA")  # alanine
end

@testset "Extract protein sequences from PDB" begin

    file(code) = joinpath(pwd(), "data", string(uppercase(code), ".pdb"))

    @testset "2VQC: Missings & selenomethionines" begin
        res = read(file("2VQC"), PDBFile)
        # seq length : 118
        # missing residues : 1-9, 80-118
        # number of modelled residues : 70
        seqs = modelled_sequences(res) # default group : All (ATOM + HETATM)
        key = (model="1", chain="A")
        seq = seqs[key]
        @test length(seq) == 70
        @test seq == "TLNSYKMAEIMYKILEKKGELTLEDILAQFEISVPSAYNIQRALKAICERHPDECEVQYKNRKTTFKWIK"
        # the first methionine is in fact a selenomethionine (MSE)
        # "HETATM   52  CA  MSE A  10"
        seq_atom = modelled_sequences(res, group="ATOM")[key]
        @test length(seq_atom) == 70 - 2 # there are two missing selenomethionines in the 
        #                                  sequence because we only selected ATOM
        #                  TLNSYK M AEI M YKI
        @test seq_atom == "TLNSYKAEIYKILEKKGELTLEDILAQFEISVPSAYNIQRALKAICERHPDECEVQYKNRKTTFKWIK"
        @test modelled_sequences(res, group="HETATM")[key] == "MM"
    end

    @testset "1AS5: NMR" begin
        res = read(file("1AS5"), PDBFile)
        seqs = modelled_sequences(res)

        @test length(seqs) == 14
        for i in 1:14
            seq = seqs[(model=string(i), chain="A")]
            @test length(seq) == 24
            @test seq == "HPPCCLYGKCRRYPGCSSASCCQR"
        end

        # test the order
        @test join(collect(k.model for k in keys(seqs))) == "1234567891011121314"

        seqs_model_1 = modelled_sequences(res, model="1")
        @test length(seqs_model_1) == 1
        @test seqs_model_1[(model="1", chain="A")] == seqs[(model="1", chain="A")]
    end

    @testset "1IGY: Insertions & chains" begin
        res = read(file("1IGY"), PDBFile)
        seqs = modelled_sequences(res)
        # chains
        @test length(seqs) == 4
        @test seqs[(model="1", chain="A")] == seqs[(model="1", chain="C")]
        @test seqs[(model="1", chain="B")] == seqs[(model="1", chain="D")]
        # insertions in chain B
        seqs_chain_B = modelled_sequences(res, chain="B")
        seq_B = seqs_chain_B[(model="1", chain="B")]
        @test seqs[(model="1", chain="B")] == seq_B
        @test occursin("LSSL", seq_B)
    end
end