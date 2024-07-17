struct PIR <: MSAFormat end

# PIR Parser
# ==========

# Description of PIR/NBRF format in :
# - https://salilab.org/modeller/9v7/manual/node445.html
# - http://iubio.bio.indiana.edu/soft/molbio/readseq/classic/src/Formats
# - http://www.bioinformatics.nl/tools/crab_pir.html
# - http://emboss.sourceforge.net/docs/themes/seqformats/NbrfFormat.html

# They could be annotations/comments after the sequence, we are ignoring that at the moment.

function _pre_readpir(io::Union{IO,AbstractString})
    IDS = String[]
    SEQS = String[]
    GS = Dict{Tuple{String,String},String}()

    # Sequence types from http://iubio.bio.indiana.edu/soft/molbio/readseq/classic/src/Formats
    nucleic_types = Set([
        "DL",  # DNA, linear
        "DC",  # DNA, circular
        "RL",  # RNA, linear
        "RC",  # RNA, circular
        "N1",  # functional RNA, other than tRNA
        "N3",
    ]) # tRNA

    nonres = r"[^A-Za-z-\.]+" # large negated character class becuase seqs are free text

    finished = false
    line_number = 0
    seq_number = 0
    for line::String in lineiterator(io)
        line_number += 1
        if startswith(line, '>')
            line_number = 1
            seq_number += 1
            finished = false
            m = match(r"^>([A-Z][A-Z0-9]);(.+)$", line) # e.g. >P1;5fd1
            if m !== nothing
                seq_type = m[1] === nothing ? "" : m[1]
                seq_id = m[2] === nothing ? "" : rstrip(m[2]::SubString)
                push!(IDS, seq_id)
                push!(GS, (seq_id, "Type") => seq_type)
                if seq_type in nucleic_types
                    @warn("Type $seq_type indicates that $seq_id is a nucleic acid.")
                end
            else
                throw(ErrorException("Identifier doesn't match the PIR format: $line"))
            end
        elseif line_number == 2 && !isempty(IDS)
            push!(GS, (IDS[seq_number], "Title") => line)
        elseif line_number == 3
            if endswith(line, '*')
                finished = true
            end
            push!(SEQS, replace(line, nonres => ""))
        elseif !finished && length(SEQS) == seq_number && seq_number != 0
            if endswith(line, '*')
                finished = true
            end
            @inbounds SEQS[seq_number] = SEQS[seq_number] * replace(line, nonres => "")
        end
    end

    (IDS, SEQS, GS)
end

function _load_sequences(
    io::Union{IO,AbstractString},
    format::Type{PIR};
    create_annotations::Bool = false,
)
    IDS, SEQS, GS = _pre_readpir(io)
    annot = Annotations()
    if create_annotations
        annot.sequences = GS
    end
    return IDS, SEQS, annot
end

# Print Pfam
# ==========

function _print_pir_seq(io::IO, seq_type, seq_id, seq_title, seq)
    println(io, ">", seq_type, ';', seq_id)
    println(io, seq_title)
    char_counter = 0
    for char in seq
        print(io, char)
        char_counter += 1
        if char_counter % 80 == 0
            println(io)
        end
    end
    println(io, '*')
end

function _get_pir_annotations(sequence_annotations, seq_id::String)
    if haskey(sequence_annotations, (seq_id, "Type"))
        seq_type = sequence_annotations[(seq_id, "Type")]
    else
        @warn("There is not sequence Type annotation for $seq_id, using XX (Unknown).")
        seq_type = "XX"
    end
    if haskey(sequence_annotations, (seq_id, "Title"))
        seq_title = sequence_annotations[(seq_id, "Title")]
    else
        @warn("There is not sequence Title annotation for $seq_id.")
        seq_title = ""
    end
    seq_type, seq_title
end

function Utils.print_file(
    io::IO,
    msa::AnnotatedMultipleSequenceAlignment,
    format::Type{PIR},
)
    sequence_annotations = getannotsequence(msa)
    seqnames = sequencenames(msa)
    for i = 1:nsequences(msa)
        seq_id = seqnames[i]
        seq_type, seq_title = _get_pir_annotations(sequence_annotations, seq_id)
        seq = stringsequence(msa, i)
        _print_pir_seq(io, seq_type, seq_id, seq_title, seq)
    end
end

function Utils.print_file(io::IO, msa::AbstractMatrix{Residue}, format::Type{PIR})
    seqnames = sequencenames(msa)
    for i = 1:nsequences(msa)
        _print_pir_seq(io, "XX", seqnames[i], "", stringsequence(msa, i))
    end
end
