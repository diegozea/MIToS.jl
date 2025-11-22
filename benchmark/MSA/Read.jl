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
    SUITE["MSA"]["read"]["FASTA.gz_annotated"] =
        @benchmarkable read_file($fasta_gz, FASTA, AnnotatedMultipleSequenceAlignment)

    global function create_gappy_fasta()
        dir = mktempdir()
        gappy_fasta = joinpath(dir, "gappy.fasta")
        open(gappy_fasta, "w") do io
            seq = repeat("A-HAGLLSAPGSCW------", 300)
            # it generates full-gap columns, both blocks and isolated gaps
            for i = 1:200
                println(io, ">seq", i)
                println(io, seq)
            end
        end
        return gappy_fasta
    end

    SUITE["MSA"]["read"]["FASTA_deletefullgaps"] = @benchmarkable read_file(
        gappy_fasta,
        FASTA,
        MultipleSequenceAlignment,
    ) setup=(gappy_fasta = create_gappy_fasta()) teardown=(
        rm(dirname(gappy_fasta); recursive = true, force = true)
    )
    SUITE["MSA"]["read"]["FASTA_deletefullgaps_mapping"] = @benchmarkable read_file(
        gappy_fasta,
        FASTA,
        AnnotatedMultipleSequenceAlignment,
        generatemapping = true,
    ) setup=(gappy_fasta = create_gappy_fasta()) teardown=(
        rm(dirname(gappy_fasta); recursive = true, force = true)
    )
end
