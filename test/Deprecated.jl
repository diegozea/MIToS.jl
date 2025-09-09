@testset "Deprecated" begin
    seq1 = res"ACDE"
    seq2 = res"ABCD"
    msa_matrix = permutedims(hcat(seq1, seq2))
    ann_msa = AnnotatedMultipleSequenceAlignment(msa_matrix)
    msa = MultipleSequenceAlignment(msa_matrix)
    pdbfile = joinpath(DATA, "short.pdb")
    pdb = read_file(pdbfile, PDBFile)

    @testset "IO" begin
        tmp1 = tempname()
        write_file(tmp1, msa_matrix, Stockholm)
        tmp2 = tempname()
        @test_deprecated write(tmp2, msa_matrix, Stockholm)
        @test read(tmp1, String) == read(tmp2, String)
        rm(tmp1)
        rm(tmp2)

        io1 = IOBuffer()
        print_file(io1, msa_matrix, Stockholm)
        str1 = String(take!(io1))
        io2 = IOBuffer()
        @test_deprecated print(io2, msa_matrix, Stockholm)
        str2 = String(take!(io2))
        @test str1 == str2

        tmp1 = tempname()
        write_file(tmp1, msa_matrix, Stockholm)
        msa_r = read_file(tmp1, Stockholm)
        @test_deprecated read(tmp1, Stockholm) == msa_r
        str_msa = read(tmp1, String)
        @test_deprecated parse(str_msa, Stockholm) == parse_file(str_msa, Stockholm)
        rm(tmp1)
    end

    @testset "MSA" begin
        pair = 1:2 .=> 1:2
        joined_new = join_msas(ann_msa, ann_msa, pair)
        @test_deprecated join(ann_msa, ann_msa, pair) == joined_new

        rng = MersenneTwister(1)
        copy1 = copy(msa_matrix)
        copy2 = copy(msa_matrix)
        expected = shuffle_msa!(rng, copy1; dims = 1, fixedgaps = true) |> getresidues
        rng = MersenneTwister(1)
        @test_deprecated shuffle!(rng, copy2, 1, true) == expected

        rng = MersenneTwister(1)
        expected2 = shuffle_msa(rng, msa_matrix; dims = 2, fixedgaps = false) |> getresidues
        rng = MersenneTwister(1)
        @test_deprecated shuffle(rng, msa_matrix, 2, false) == expected2

        @test_deprecated convert(MultipleSequenceAlignment, ann_msa) ==
                         MultipleSequenceAlignment(ann_msa)
        seq_ann = AnnotatedAlignedSequence(msa_matrix[1:1, :])
        @test_deprecated convert(AlignedSequence, seq_ann) == AlignedSequence(seq_ann)
        @test_deprecated convert(AnnotatedMultipleSequenceAlignment, msa) ==
                         AnnotatedMultipleSequenceAlignment(msa)
        seq_plain = AlignedSequence(msa_matrix[1:1, :])
        @test_deprecated convert(AnnotatedAlignedSequence, seq_plain) ==
                         AnnotatedAlignedSequence(seq_plain)

        @test_deprecated transpose(seq_plain) == permutedims(seq_plain)
    end

    @testset "Statistics" begin
        s = res"ARNDCQEGHILKMFPSTWYV-"
        Ps = probabilities(s)
        @test_deprecated entropy(Ps) == shannon_entropy(Ps)
        @test_deprecated entropy(Ps, 2) == shannon_entropy(Ps; base = 2)
        @test_deprecated marginal_entropy(Ps, 1) == marginal_entropy(Ps; margin = 1)
        @test_deprecated marginal_entropy(Ps, 1, 2) ==
                         marginal_entropy(Ps; margin = 1, base = 2)
        @test_deprecated kullback_leibler(Ps, Ps, 2) ==
                         kullback_leibler(Ps; background = Ps, base = 2)
        @test_deprecated kullback_leibler(Ps, Ps) == kullback_leibler(Ps; background = Ps)
        @test_deprecated kullback_leibler(Ps, 2) == kullback_leibler(Ps; base = 2)
        @test_deprecated mutual_information(probabilities(s, s), 2) ==
                         mutual_information(probabilities(s, s); base = 2)
    end

    @testset "Frequencies" begin
        table = ContingencyTable(Float64, Val{1}, UngappedAlphabet())
        res_counts = (@test_deprecated Counts{Float64,1,UngappedAlphabet}(table))
        @test isa(res_counts, Frequencies{Float64,1,UngappedAlphabet})

        aln = read_file(joinpath(DATA, "Gaoetal2011.fasta"), FASTA)
        t_old = Frequencies(ContingencyTable(Float64, Val{2}, UngappedAlphabet()))
        t_new = Frequencies(ContingencyTable(Float64, Val{2}, UngappedAlphabet()))
        res_new =
            mapcolpairfreq!(normalized_mutual_information, aln, t_new; usediagonal = false)
        res_old = @test_deprecated mapcolpairfreq!(
            normalized_mutual_information,
            aln,
            t_old,
            Val{false},
        )
        @test getarray(res_old) == getarray(res_new)

        t_old2 = Frequencies(ContingencyTable(Float64, Val{2}, UngappedAlphabet()))
        t_new2 = Frequencies(ContingencyTable(Float64, Val{2}, UngappedAlphabet()))
        res_new2 = mapseqpairfreq!(
            normalized_mutual_information,
            permutedims(aln),
            t_new2;
            usediagonal = false,
        )
        res_old2 = @test_deprecated mapseqpairfreq!(
            normalized_mutual_information,
            permutedims(aln),
            t_old2,
            Val{false},
        )
        @test getarray(res_old2) == getarray(res_new2)

        seqs = (seq1,)
        tab1 = ContingencyTable(Float64, Val{1}, UngappedAlphabet())
        tab2 = ContingencyTable(Float64, Val{1}, UngappedAlphabet())
        frequencies!(tab1, seqs...)
        @test_deprecated count!(tab2, NoClustering(), NoPseudocount(), seqs...)
        @test tab1 == tab2

        @test_deprecated count(seq1; alphabet = UngappedAlphabet()) ==
                         frequencies(seq1; alphabet = UngappedAlphabet())

        pt1 = ContingencyTable(Float64, Val{1}, UngappedAlphabet())
        pt2 = ContingencyTable(Float64, Val{1}, UngappedAlphabet())
        probabilities!(
            pt1,
            seqs...;
            weights = NoClustering(),
            pseudocounts = NoPseudocount(),
            pseudofrequencies = NoPseudofrequencies(),
        )
        @test_deprecated probabilities!(
            pt2,
            NoClustering(),
            NoPseudocount(),
            NoPseudofrequencies(),
            seqs...,
        )
        @test pt1 == pt2

        msa_file = joinpath(DATA, "Gaoetal2011.fasta")
        res_new_buslje = buslje09(
            read_file(msa_file, FASTA);
            samples = 0,
            clustering = false,
            apc = false,
        )
        @test_deprecated buslje09(
            msa_file,
            FASTA;
            samples = 0,
            clustering = false,
            apc = false,
        ) == res_new_buslje
        res_new_blmi = BLMI(read_file(msa_file, FASTA); samples = 0, apc = false)
        @test_deprecated BLMI(msa_file, FASTA; samples = 0, apc = false) == res_new_blmi
    end

    @testset "PDB" begin
        res = pdb[1]
        @test_deprecated isresidue(res.id, "1", "A", "ATOM", "1") == isresidue(
            res.id;
            model = "1",
            chain = "A",
            group = "ATOM",
            residue = "1",
        )
        @test_deprecated isresidue(res, "1", "A", "ATOM", "1") == isresidue(
            res;
            model = "1",
            chain = "A",
            group = "ATOM",
            residue = "1",
        )
        @test_deprecated residues(pdb, "1", "A", "ATOM", "1") == select_residues(
            pdb;
            model = "1",
            chain = "A",
            group = "ATOM",
            residue = "1",
        )
        @test_deprecated residuesdict(pdb, "1", "A", "ATOM", "1") == residuesdict(
            pdb;
            model = "1",
            chain = "A",
            group = "ATOM",
            residue = "1",
        )
        @test_deprecated eval(
            :(@residues $pdb model "1" chain "A" group "ATOM" residue "1"),
        ) == select_residues(
            pdb;
            model = "1",
            chain = "A",
            group = "ATOM",
            residue = "1",
        )
        @test_deprecated eval(
            :(@residuesdict $pdb model "1" chain "A" group "ATOM" residue "1"),
        ) == residuesdict(
            pdb;
            model = "1",
            chain = "A",
            group = "ATOM",
            residue = "1",
        )
        @test_deprecated atoms(pdb, "1", "A", "ATOM", "1", "CA") == select_atoms(
            pdb;
            model = "1",
            chain = "A",
            group = "ATOM",
            residue = "1",
            atom = "CA",
        )
        @test_deprecated eval(
            :(@atoms $pdb model "1" chain "A" group "ATOM" residue "1" atom "CA"),
        ) == select_atoms(
            pdb;
            model = "1",
            chain = "A",
            group = "ATOM",
            residue = "1",
            atom = "CA",
        )
    end

    @testset "Utilities" begin
        vec = [1, 2, 3]
        @test_deprecated deleteitems!(vec, [2])
        @test vec == [1, 3]
    end
end
