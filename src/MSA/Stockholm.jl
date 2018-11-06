struct Stockholm <: FileFormat end

@inline function _fill_with_sequence_line!(IDS, SEQS, line)
    if !startswith(line,'#') && !startswith(line,"//")
        words = get_n_words(line, 2)
        @inbounds id = words[1]
        if id in IDS
            # It's useful when sequences are split into several lines
            # It can be a problem with duplicated IDs
            i = something(findfirst(isequal(id), IDS), 0)
            SEQS[i] = SEQS[i] * words[2]
        else
            push!(IDS, id)
            push!(SEQS, words[2])
        end
    end
end

function _fill_with_line!(IDS, SEQS, GF, GS, GC, GR, line)
    if startswith(line,"#=GF")
        words = get_n_words(line, 3)
        id = words[2]
        if id in keys(GF)
            GF[ id ] = GF[ id ] * "\n" * words[3]
        else
            GF[ id ] = words[3]
        end
    elseif startswith(line,"#=GS")
        words = get_n_words(line, 4)
        idtuple = (words[2],words[3])
        if idtuple in keys(GS)
            GS[ idtuple ] = GS[ idtuple ] * "\n" * words[4]
        else
            GS[ idtuple ] = words[4]
        end
    elseif startswith(line,"#=GC")
        words = get_n_words(line, 3)
        GC[words[2]] = words[3]
    elseif startswith(line,"#=GR")
        words = get_n_words(line, 4)
        GR[(words[2],words[3])] = words[4]
    else
        _fill_with_sequence_line!(IDS, SEQS, line)
    end
end

function _pre_readstockholm(io::Union{IO, AbstractString})
    IDS  = OrderedSet{String}()
    SEQS = String[]
    GF = OrderedDict{String,String}()
    GC = Dict{String,String}()
    GS = Dict{Tuple{String,String},String}()
    GR = Dict{Tuple{String,String},String}()

    @inbounds for line::String in lineiterator(io)
        _fill_with_line!(IDS, SEQS, GF, GS, GC, GR, line)
        if startswith(line, "//")
           break
        end
    end

    _check_seq_len(IDS, SEQS)

    GF = sizehint!(GF, length(GF))
    GC = sizehint!(GC, length(GC))
    GS = sizehint!(GS, length(GS))
    GR = sizehint!(GR, length(GR))
    (IDS, SEQS, GF, GS, GC, GR)
end

function _pre_readstockholm_sequences(io::Union{IO, AbstractString})
    IDS  = OrderedSet{String}()
    SEQS = String[]
    @inbounds for line::String in lineiterator(io)
        _fill_with_sequence_line!(IDS, SEQS, line)
        if startswith(line, "//")
           break
        end
    end
    _check_seq_len(IDS, SEQS)
    (IDS, SEQS)
end

function Base.parse(io::Union{IO, AbstractString},
                   format::Type{Stockholm},
                   output::Type{AnnotatedMultipleSequenceAlignment};
                   generatemapping::Bool=false,
                   useidcoordinates::Bool=false,
                   deletefullgaps::Bool=true,
                   keepinserts::Bool=false)::AnnotatedMultipleSequenceAlignment
    IDS, SEQS, GF, GS, GC, GR = _pre_readstockholm(io)
    annot = Annotations(GF, GS, GC, GR)
    _generate_annotated_msa(annot, collect(IDS), SEQS, keepinserts,
                            generatemapping, useidcoordinates, deletefullgaps)
end

function Base.parse(io::Union{IO, AbstractString},
                   format::Type{Stockholm},
                   output::Type{NamedResidueMatrix{Array{Residue,2}}};
                   deletefullgaps::Bool=true)::NamedResidueMatrix{Array{Residue,2}}
    IDS, SEQS = _pre_readstockholm_sequences(io)
    msa = _generate_named_array(SEQS, IDS)
    if deletefullgaps
        return deletefullgapcolumns(msa)
    end
    msa
end

function Base.parse(io::Union{IO, AbstractString},
                   format::Type{Stockholm},
                   output::Type{MultipleSequenceAlignment};
                   deletefullgaps::Bool=true)::MultipleSequenceAlignment
    msa = parse(io,
                format,
                NamedResidueMatrix{Array{Residue,2}},
                deletefullgaps=deletefullgaps)
    MultipleSequenceAlignment(msa)
end

function Base.parse(io::Union{IO,AbstractString},
                   format::Type{Stockholm},
                   output::Type{Matrix{Residue}};
                   deletefullgaps::Bool=true)::Matrix{Residue}
    IDS, SEQS = _pre_readstockholm_sequences(io)
    _strings_to_matrix_residue_unsafe(SEQS, deletefullgaps)
end

function Base.parse(io, format::Type{Stockholm};
                    generatemapping::Bool=false,
                    useidcoordinates::Bool=false,
                    deletefullgaps::Bool=true,
                    keepinserts::Bool=false)::AnnotatedMultipleSequenceAlignment
    parse(io, Stockholm, AnnotatedMultipleSequenceAlignment,
          generatemapping=generatemapping,
          useidcoordinates=useidcoordinates,
          deletefullgaps=deletefullgaps,
          keepinserts=keepinserts)
end

# Print Pfam
# ==========

function _to_sequence_dict(annotation::Dict{Tuple{String,String},String})
    seq_dict = Dict{String,Vector{String}}()
    for (key, value) in annotation
        seq_id = key[1]
        if haskey(seq_dict, seq_id)
            push!(seq_dict[seq_id], string(seq_id, '\t', key[2], '\t', value))
        else
            seq_dict[seq_id] = [ string(seq_id, '\t', key[2], '\t', value) ]
        end
    end
    sizehint!(seq_dict, length(seq_dict))
end

function Base.print(io::IO, msa::AnnotatedMultipleSequenceAlignment,
                    format::Type{Stockholm})
    _printfileannotations(io, msa.annotations)
    _printsequencesannotations(io, msa.annotations)
    res_annotations = _to_sequence_dict(msa.annotations.residues)
    seqnames = sequencenames(msa)
    for i in 1:nsequences(msa)
        id = seqnames[i]
        seq = stringsequence(msa, i)
        println(io, id, "\t\t\t", seq)
        if id in keys(res_annotations)
            for line in res_annotations[id]
                println(io, "#=GR ", line)
            end
        end
    end
    _printcolumnsannotations(io, msa.annotations)
    println(io, "//")
end

function Base.print(io::IO, msa::AbstractMatrix{Residue}, format::Type{Stockholm})
    seqnames = sequencenames(msa)
    for i in 1:nsequences(msa)
        println(io, seqnames[i], "\t\t\t", stringsequence(msa, i))
    end
    println(io, "//")
end

Base.print(msa::AnnotatedMultipleSequenceAlignment) = print(stdout, msa, Stockholm)
