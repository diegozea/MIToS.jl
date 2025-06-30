let
    msa_file = joinpath(@__DIR__, "..", "..", "test", "data", "simple.fasta")
    msa_a = read_file(
        msa_file,
        FASTA,
        AnnotatedMultipleSequenceAlignment,
        generatemapping = true,
    )
    msa_b = read_file(
        msa_file,
        FASTA,
        AnnotatedMultipleSequenceAlignment,
        generatemapping = true,
    )
    SUITE["MSA"]["Base.vcat"]["annotated"] = @benchmarkable vcat($msa_a, $msa_b)
    msa_u_a = MultipleSequenceAlignment(msa_a)
    msa_u_b = MultipleSequenceAlignment(msa_b)
    SUITE["MSA"]["Base.vcat"]["unannotated"] = @benchmarkable vcat($msa_u_a, $msa_u_b)
end
