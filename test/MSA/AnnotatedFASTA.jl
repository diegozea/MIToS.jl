@testset "AnnotatedFASTASequences" begin

    @testset "_parse_aff_header with README sample header" begin
        # https://github.com/NawarMalhis/AFF/blob/f0e7033bc4d08b255cd883c99523d9379a9fbb16/README.md
        header_lines = [
            "# Data Name: A sample annotated FASTA format file",
            "#",
            "# Optional top comments line 1",
            "# Optional top comments line 2",
            "#",
            "# Sequences: 34",
            "#",
            "# Format:",
            "# >accession|DisProt=DisProt_ID|UniProt=UniProt_ID|OX=OX_ID|Pfam=Pfam_ID|IDEAL=IDEAL_ID|PDB=PDB_ID|AlphaFoldDB=AlphaFoldDB_ID|ELM=ELM_ID",
            "# Amino acid sequence",
            "# IDR  Protein disordered region",
            "# DtO  Disordered to ordered transition",
            "# Linker  Linker regions",
            "# binding  IDR binding in general",
            "# binding_protein  IDR-binding to proteins",
            "# binding_nucleic  IDR-binding to nucleic",
            "# binding_ion  IDR-binding to ions",
            "# binding_lipid  IDR-binding to lipids",
            "# binding_SM  IDR binding to small molecules",
            "#",
            "# ID Counts:",
            "#   ID   Seq#   Total#  Unique#",
            "#   DisProt  34  34   34",
            "#   UniProt  34  137  137",
            "#   OX  34  141  99",
            "#   Pfam  33  58  56",
            "#   IDEAL  8  8  8",
            "#   PDB  29  509  509",
            "#   AlphaFoldDB  31  95  95",
            "#   ELM  4  4  4",
            "#",
            "# Tags Counts:",
            "#   tag  Seq# Seg#  '0' '1'  '-'",
            "#   IDR  34  44  8,016  4,247  11",
            "#   DtO  12  16  1,803  775  0",
            "#   Linker  6  6  1,381  93  0",
            "#   binding  18  19  3,588  2,534  0",
            "#   binding_protein  13  14  2,674  2,138  0",
            "#   binding_nucleic  3  3  44  240  0",
            "#   binding_ion  0  0  0  0  0",
            "#   binding_lipid  3  3  551  361  0",
            "#   binding_SM  1  1  319  9  0",
            "#",
            "# Optional bottom comments line 1",
            "# Optional bottom comments line 2",
            "#",
        ]

        file_annotations, format_tags = MIToS.MSA._parse_aff_header(header_lines)

        @test file_annotations["Data Name"] == "A sample annotated FASTA format file"
        @test file_annotations["Sequences"] == "34"
        @test file_annotations["Comment_1"] == "Optional top comments line 1"
        @test file_annotations["Comment_2"] == "Optional top comments line 2"
        @test file_annotations["Comment_3"] == "Optional bottom comments line 1"
        @test file_annotations["Comment_4"] == "Optional bottom comments line 2"
        @test haskey(file_annotations, "Format")
        @test haskey(file_annotations, "ID Counts")
        @test haskey(file_annotations, "Tags Counts")
        @test occursin(">accession|DisProt=DisProt_ID", file_annotations["Format"])
        @test occursin(
            "binding_SM  IDR binding to small molecules",
            file_annotations["Format"],
        )
        @test occursin("DisProt  34  34   34", file_annotations["ID Counts"])
        @test occursin("ELM  4  4  4", file_annotations["ID Counts"])
        @test occursin("tag  Seq# Seg#  '0' '1'  '-'", file_annotations["Tags Counts"])
        @test occursin("binding_lipid  3  3  551  361  0", file_annotations["Tags Counts"])
        @test format_tags == [
            "IDR",
            "DtO",
            "Linker",
            "binding",
            "binding_protein",
            "binding_nucleic",
            "binding_ion",
            "binding_lipid",
            "binding_SM",
        ]
    end

    @testset "_parse_aff_header section boundary without blank line" begin
        header_lines = [
            "# Data Name: Minimal",
            "# Format:",
            "# >accession|FASTA=identifier",
            "# Amino acid sequence",
            "# IDR  Protein disordered region",
            "# ID Counts:",
            "#   ID   Seq#   Total#  Unique#",
            "#   FASTA  1  1  1",
        ]

        file_annotations, format_tags = MIToS.MSA._parse_aff_header(header_lines)

        @test file_annotations["Data Name"] == "Minimal"
        @test haskey(file_annotations, "Format")
        @test haskey(file_annotations, "ID Counts")
        @test occursin("IDR  Protein disordered region", file_annotations["Format"])
        @test occursin("FASTA  1  1  1", file_annotations["ID Counts"])
        @test format_tags == ["IDR"]
    end

    @testset "format_has_sequence_line karg in read_file" begin
        with_descriptor = """
# Format:
# >accession|FASTA=identifier
# RNA seq
# IDR  Protein disordered region
>SEQ1
AC
01
"""
        mktemp() do path, io
            write(io, with_descriptor)
            close(io)
            parsed = read_file(path, AnnotatedFASTASequences)
            @test getannotresidue(parsed[1], "IDR") == "01"
        end

        # If Format omits descriptor line, set `format_has_sequence_line=false`.
        without_descriptor = """
# Format:
# >accession|FASTA=identifier
# IDR  Protein disordered region
>SEQ2
AC
01
"""
        mktemp() do path, io
            write(io, without_descriptor)
            close(io)
            parsed =
                read_file(path, AnnotatedFASTASequences; format_has_sequence_line = false)
            @test getannotresidue(parsed[1], "IDR") == "01"
        end
    end

    @testset "Format tags count mismatch throws informative error" begin
        mismatched = """
# Format:
# >accession|FASTA=identifier
# Amino acid sequence
# IDR  Protein disordered region
# binding  Binding region
>SEQ_MISMATCH
AC
01
"""

        err = try
            parse_file(IOBuffer(mismatched), AnnotatedFASTASequences)
            nothing
        catch e
            e
        end

        @test err isa ErrorException
        @test occursin("The sequence SEQ_MISMATCH has 1 annotation lines", err.msg)
        @test occursin("the header Format section defines 2 features", err.msg)
    end

    @testset "Missing sequence and annotation lines throw informative errors" begin
        no_sequence_lines = """
>NO_SEQUENCE
"""

        err_no_seq = try
            parse_file(IOBuffer(no_sequence_lines), AnnotatedFASTASequences)
            nothing
        catch e
            e
        end

        @test err_no_seq isa ErrorException
        @test occursin("There are no sequence lines for NO_SEQUENCE", err_no_seq.msg)

        no_annotation_lines = """
>NO_ANNOT
ACDE
"""

        err_no_annot = try
            parse_file(IOBuffer(no_annotation_lines), AnnotatedFASTASequences)
            nothing
        catch e
            e
        end

        @test err_no_annot isa ErrorException
        @test occursin("There are no annotation lines for NO_ANNOT", err_no_annot.msg)
    end

    @testset "Annotation length mismatch throws informative error" begin
        length_mismatch = """
>LEN_MISMATCH
ACDE
01
"""

        err = try
            parse_file(IOBuffer(length_mismatch), AnnotatedFASTASequences)
            nothing
        catch e
            e
        end

        @test err isa ErrorException
        @test occursin("The annotation line 1 for LEN_MISMATCH has 2 characters", err.msg)
        @test occursin("4 are expected", err.msg)
    end

    @testset "Duplicate sequence names keep OriginalSeqName annotation" begin
        duplicated_ids = """
>dup
AC
01
>dup
GT
10
"""

        data = parse_file(IOBuffer(duplicated_ids), AnnotatedFASTASequences)
        @test length(data) == 2
        @test sequence_id(data[1]) == "dup"
        @test sequence_id(data[2]) == "dup(1)"
        @test getannotsequence(data[2], "OriginalSeqName") == "dup"
    end

    @testset "AFF" begin
        aff = read_file(joinpath(DATA, "aff_example.fasta"), AnnotatedFASTASequences)

        @test isa(aff, Vector{AnnotatedSequence})
        @test length(aff) == 2

        seq1 = aff[1]
        @test sequence_id(seq1) == "AFF_SEQ_1|FASTA=AFF_SEQ_1"
        @test join(seq1) == "ACDEFG"
        @test getannotfile(seq1, "Data Name") == "Example AFF file"
        @test getannotfile(seq1, "Sequences") == "2"
        @test getannotresidue(seq1, "IDR") == "010011"
        @test getannotresidue(seq1, "binding") == "111000"

        seq2 = aff[2]
        @test sequence_id(seq2) == "AFF_SEQ_2|FASTA=AFF_SEQ_2"
        @test join(seq2) == "GGTTAA"
        @test getannotresidue(seq2, "IDR") == "111000"
        @test getannotresidue(seq2, "binding") == "000111"
    end

    @testset "CAID-like (no header, one annotation line)" begin
        caid = read_file(joinpath(DATA, "caid_modified.fasta"), AnnotatedFASTASequences)

        @test isa(caid, Vector{AnnotatedSequence})
        @test length(caid) == 2

        seq1 = caid[1]
        @test sequence_id(seq1) == "CAID_SEQ_1"
        @test join(seq1) == "ACDEFGH"
        @test isempty(getannotfile(seq1))
        @test isempty(getannotsequence(seq1))
        @test getannotresidue(seq1, "Feature") == "0100110"

        seq2 = caid[2]
        @test sequence_id(seq2) == "CAID_SEQ_2"
        @test join(seq2) == "MNPQRS"
        @test isempty(getannotsequence(seq2))
        @test getannotresidue(seq2, "Feature") == "11-000"
    end

    @testset "No header, multiple annotation lines" begin
        noheader =
            read_file(joinpath(DATA, "aff_headerless.fasta"), AnnotatedFASTASequences)

        @test isa(noheader, Vector{AnnotatedSequence})
        @test length(noheader) == 2

        seq1 = noheader[1]
        @test sequence_id(seq1) == "MULTI_SEQ_1"
        @test join(seq1) == "ACDEFG"
        @test getannotresidue(seq1, "Feature_1") == "010011"
        @test getannotresidue(seq1, "Feature_2") == "111000"

        seq2 = noheader[2]
        @test sequence_id(seq2) == "MULTI_SEQ_2"
        @test join(seq2) == "MNPQRS"
        @test getannotresidue(seq2, "Feature_1") == "00-000"
        @test getannotresidue(seq2, "Feature_2") == "101010"
    end

    @testset "CAID-like annotations preserve trailing spaces" begin
        labels = "01" * repeat(" ", 4) # user-defined unknown labels as trailing spaces
        content = ">SPACE_SEQ\nACDEFG\n" * labels * "\n"

        mktemp() do path, io
            write(io, content)
            close(io)

            data = read_file(path, AnnotatedFASTASequences)
            @test length(data) == 1
            @test join(data[1]) == "ACDEFG"
            @test getannotresidue(data[1], "Feature") == labels
        end
    end
end
