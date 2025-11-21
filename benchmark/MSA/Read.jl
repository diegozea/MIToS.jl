let
    data_dir = joinpath(@__DIR__, "..", "..", "test", "data")
    fasta_gz = joinpath(data_dir, "PF09645_full.fasta.gz")
    sth = joinpath(data_dir, "PF09645_full.stockholm")
    gappy_fasta = begin
        dir = mktempdir()
        atexit(() -> rm(dir; force = true, recursive = true))
        path = joinpath(dir, "gappy_read.fasta")
        open(path, "w") do io
            seq = repeat("A", 300) * repeat("-", 100)
            for i in 1:200
                println(io, ">seq", i)
                println(io, seq)
            end
        end
        path
    end

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
    SUITE["MSA"]["read"]["FASTA_deletefullgaps"] = @benchmarkable read_file(
        $gappy_fasta,
        FASTA,
        Matrix{Residue},
    )
end
