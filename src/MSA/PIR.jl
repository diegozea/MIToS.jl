struct PIR <: FileFormat end

# PIR Parser
# ==========

# Description of PIR/NBRF format in :
# - https://salilab.org/modeller/9v7/manual/node445.html
# - http://iubio.bio.indiana.edu/soft/molbio/readseq/classic/src/Formats
# - http://www.bioinformatics.nl/tools/crab_pir.html
# - http://emboss.sourceforge.net/docs/themes/seqformats/NbrfFormat.html

# They could be annotations/comments after the sequence, we are ignoring that at the moment.

function _pre_readpir(io::Union{IO, AbstractString})
    IDS = String[]
    SEQS = String[]
    GS = Dict{Tuple{String,String},String}()

    # Sequence types from http://iubio.bio.indiana.edu/soft/molbio/readseq/classic/src/Formats
    nucleic_types = Set(["DL",  # DNA, linear
                         "DC",  # DNA, circular
                         "RL",  # RNA, linear
                         "RC",  # RNA, circular
                         "N1",  # functional RNA, other than tRNA
                         "N3"]) # tRNA

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
            m = match(r"^>([A-Z][A-Z0-9]);(\S+)", line) # e.g. >P1;5fd1
            if m !== nothing
                seq_type = m[1]
                seq_id   = m[2]
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

function Base.parse(io::Union{IO, AbstractString},
                   format::Type{PIR},
                   output::Type{AnnotatedMultipleSequenceAlignment};
                   generatemapping::Bool=false,
                   useidcoordinates::Bool=false,
                   deletefullgaps::Bool=true,
                   keepinserts::Bool=false)::AnnotatedMultipleSequenceAlignment
    IDS, SEQS, GS = _pre_readpir(io)
    annot = Annotations()
    annot.sequences = GS
    _generate_annotated_msa(annot, IDS, SEQS, keepinserts, generatemapping,
                            useidcoordinates, deletefullgaps)
end

function Base.parse(io::Union{IO, AbstractString},
                   format::Type{PIR},
                   output::Type{NamedResidueMatrix{Array{Residue,2}}};
                   deletefullgaps::Bool=true)::NamedResidueMatrix{Array{Residue,2}}
    IDS, SEQS, _ = _pre_readpir(io)
    msa = _generate_named_array(SEQS, IDS)
    if deletefullgaps
        return deletefullgapcolumns(msa)
    end
    msa
end

function Base.parse(io::Union{IO, AbstractString},
                   format::Type{PIR},
                   output::Type{MultipleSequenceAlignment};
                   deletefullgaps::Bool=true)::MultipleSequenceAlignment
    msa = parse(io,
                format,
                NamedResidueMatrix{Array{Residue,2}},
                deletefullgaps=deletefullgaps)
    MultipleSequenceAlignment(msa)
end

function Base.parse(io::Union{IO,AbstractString},
                   format::Type{PIR},
                   output::Type{Matrix{Residue}};
                   deletefullgaps::Bool=true)::Matrix{Residue}
    IDS, SEQS, _ = _pre_readpir(io)
    _strings_to_matrix_residue_unsafe(SEQS, deletefullgaps)
end

function Base.parse(io, format::Type{PIR};
                    generatemapping::Bool=false,
                    useidcoordinates::Bool=false,
                    deletefullgaps::Bool=true,
                    keepinserts::Bool=false)::AnnotatedMultipleSequenceAlignment
    parse(io, PIR, AnnotatedMultipleSequenceAlignment,
          generatemapping=generatemapping,
          useidcoordinates=useidcoordinates,
          deletefullgaps=deletefullgaps,
          keepinserts=keepinserts)
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

function Base.print(io::IO, msa::AnnotatedMultipleSequenceAlignment,
                    format::Type{PIR})
    seqann = getannotsequence(msa)
    seqnames = sequencenames(msa)
    for i in 1:nsequences(msa)
        seq_id = seqnames[i]
        if haskey(seqann, (seq_id, "Type"))
            seq_type = seqann[(seq_id, "Type")]
        else
            @warn("There is not sequence Type annotation for $seq_id, using XX (Unknown).")
            seq_type = "XX"
        end
        if haskey(seqann, (seq_id, "Title"))
            seq_title = seqann[(seq_id, "Title")]
        else
            @warn("There is not sequence Title annotation for $seq_id.")
            seq_title = ""
        end
        seq = stringsequence(msa, i)
        _print_pir_seq(io, seq_type, seq_id, seq_title, seq)
    end
end

function Base.print(io::IO, msa::AbstractMatrix{Residue}, format::Type{PIR})
    seqnames = sequencenames(msa)
    for i in 1:nsequences(msa)
        _print_pir_seq(io, "XX", seqnames[i], "", stringsequence(msa, i))
    end
end
