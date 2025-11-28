using CodecZlib

let
    data_dir = joinpath(@__DIR__, "..", "..", "test", "data")
    fasta_gz = joinpath(data_dir, "PF09645_full.fasta.gz")
    sth = joinpath(data_dir, "PF09645_full.stockholm")
    clustal = joinpath(data_dir, "PF09645.aln")
    clustal_num = joinpath(data_dir, "PF09645.aln-num")

    global function create_stockholm_gz(stockholm_path)
        dir = mktempdir()
        stockholm_gz = joinpath(dir, "PF09645_full.stockholm.gz")
        open(stockholm_path) do src
            open(stockholm_gz, "w") do dest
                gz = CodecZlib.GzipCompressorStream(dest)
                write(gz, read(src))
                close(gz)
            end
        end
        return stockholm_gz
    end

    global function create_uncompressed_fasta(fasta_gz_path)
        dir = mktempdir()
        fasta_path = joinpath(dir, replace(basename(fasta_gz_path), ".gz" => ""))
        open(fasta_gz_path) do src
            open(fasta_path, "w") do dest
                stream = CodecZlib.GzipDecompressorStream(src)
                write(dest, read(stream))
                close(stream)
            end
        end
        return fasta_path
    end

    SUITE["MSA"]["read"]["Stockholm"] =
        @benchmarkable read_file($sth, Stockholm, MultipleSequenceAlignment)
    SUITE["MSA"]["read"]["Stockholm.gz"] = @benchmarkable read_file(
        stockholm_gz,
        Stockholm,
        MultipleSequenceAlignment,
    ) setup=(stockholm_gz = create_stockholm_gz($sth)) teardown=(rm(
        dirname(stockholm_gz);
        recursive = true,
        force = true,
    ))
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
    SUITE["MSA"]["read"]["FASTA"] = @benchmarkable read_file(
        fasta_path,
        FASTA,
        MultipleSequenceAlignment,
    ) setup=(fasta_path = create_uncompressed_fasta($fasta_gz)) teardown=(rm(
        dirname(fasta_path);
        recursive = true,
        force = true,
    ))
    SUITE["MSA"]["read"]["FASTA.gz_annotated"] =
        @benchmarkable read_file($fasta_gz, FASTA, AnnotatedMultipleSequenceAlignment)
    SUITE["MSA"]["read"]["Clustal"] =
        @benchmarkable read_file($clustal, Clustal, MultipleSequenceAlignment)
    SUITE["MSA"]["read"]["Clustal_num"] =
        @benchmarkable read_file($clustal_num, Clustal, MultipleSequenceAlignment)

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

    SUITE["MSA"]["read"]["FASTA_deletefullgaps"] =
        @benchmarkable read_file(gappy_fasta, FASTA, MultipleSequenceAlignment) setup=(
            gappy_fasta = create_gappy_fasta()
        ) teardown=(rm(dirname(gappy_fasta); recursive = true, force = true))
    SUITE["MSA"]["read"]["FASTA_deletefullgaps_mapping"] = @benchmarkable read_file(
        gappy_fasta,
        FASTA,
        AnnotatedMultipleSequenceAlignment,
        generatemapping = true,
    ) setup=(gappy_fasta = create_gappy_fasta()) teardown=(rm(
        dirname(gappy_fasta);
        recursive = true,
        force = true,
    ))
end
