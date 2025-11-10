let
    data_dir = joinpath(@__DIR__, "..", "..", "test", "data")
    fasta_gz = joinpath(data_dir, "PF09645_full.fasta.gz")
    sth = joinpath(data_dir, "PF09645_full.stockholm")

    SUITE["MSA"]["read"]["Stockholm"] =
        @benchmarkable read_file($sth, Stockholm, MultipleSequenceAlignment)
    SUITE["MSA"]["read"]["Stockholm_annotated"] =
        @benchmarkable read_file($sth, Stockholm, AnnotatedMultipleSequenceAlignment)
    SUITE["MSA"]["read"]["Stockholm_mapping"] = @benchmarkable read_file(
        $sth,
        Stockholm,
        AnnotatedMultipleSequenceAlignment,
        generatemapping = true,
    )
    SUITE["MSA"]["read"]["Stockholm_mapping_coords"] = @benchmarkable read_file(
        $sth,
        Stockholm,
        AnnotatedMultipleSequenceAlignment,
        generatemapping = true,
        useidcoordinates = true,
    )
    SUITE["MSA"]["read"]["FASTA.gz"] =
        @benchmarkable read_file($fasta_gz, FASTA, MultipleSequenceAlignment)
end
