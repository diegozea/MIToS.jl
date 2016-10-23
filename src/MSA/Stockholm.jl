immutable Stockholm <: Format end

@inline function _fill_with_sequence_line!(line, IDS, SEQS)
    if !startswith(line,'#') && !startswith(line,"//")
        words = get_n_words(line, 2)
        id = words[1]
        i = findfirst(IDS, id)
        if i == 0
            push!(IDS, id)
            push!(SEQS,words[2])
        else
            # It's useful when sequences are split into several lines
            # It can be a problem with duplicated IDs
            SEQS[i] = SEQS[i] * words[2]
        end
    end
end

function _fill_with_line!(line, IDS, SEQS, GF, GS, GC, GR)
    if startswith(line,"#=GF")
        words = get_n_words(line,3)
        id = words[2]
        if id in keys(GF)
            GF[ id ] = GF[ id ] * "\n" * words[3]
        else
            GF[ id ] = words[3]
        end
    elseif startswith(line,"#=GS")
        words = get_n_words(line,4)
        idtuple = (words[2],words[3])
        if idtuple in keys(GS)
            GS[ idtuple ] = GS[ idtuple ] * "\n" * words[4]
        else
            GS[ idtuple ] = words[4]
        end
    elseif startswith(line,"#=GC")
        words = get_n_words(line,3)
        GC[words[2]] = words[3]
    elseif startswith(line,"#=GR")
        words = get_n_words(line,4)
        GR[(words[2],words[3])] = words[4]
    else
        _fill_with_sequence_line!(line, IDS, SEQS)
    end
end

function _check_seq_len(IDS,SEQS)
    N = length(SEQS)
    if N != 0 & N > 1
        first_length = length(SEQS[1])
        for i in 1:N
            len = length(SEQS[i])
            if len != first_length
                if rem(len,first_length) != 0
                    throw(ErrorException("""
                        The sequence $(IDS[i]) has $len residues.
                        $first_length residues are expected.
                        """))
                else
                    n = div(len,first_length)
                    throw(ErrorException("""
                        The sequence $(IDS[i]) has $len residues.
                        $first_length residues are expected.
                        Please check if there are $n sequences with the same name.
                        """))
                end
            end
        end
    end
end

function _pre_readstockholm(io::Union{IO, AbstractString})
    IDS  = String[]
    SEQS = String[]
    GF = OrderedDict{String,String}()
    GC = Dict{String,String}()
    GS = Dict{Tuple{String,String},String}()
    GR = Dict{Tuple{String,String},String}()

    @inbounds for line in eachline(io)
        if length(line) >= 4
            _fill_with_line!(line, IDS, SEQS, GF, GS, GC, GR)
            if startswith(line,"//")
               break
            end
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
    IDS  = String[]
    SEQS = String[]
    @inbounds for line in eachline(io)
        if length(line) >= 4
            _fill_with_sequence_line!(line, IDS, SEQS)
            if startswith(line,"//")
               break
            end
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
                   keepinserts::Bool=false)
    IDS, SEQS, GF, GS, GC, GR = _pre_readstockholm(io)
    annot = Annotations(GF, GS, GC, GR)
    if keepinserts
        _keepinserts!(SEQS, annot)
    end
    if generatemapping
        if useidcoordinates && hascoordinates(IDS[1])
            MSA, MAP = _to_msa_mapping(SEQS, IDS)
        else
            MSA, MAP = _to_msa_mapping(SEQS)
            setnames!(MSA, IDS, 1)
        end
        setannotfile!(annot, "NCol", string(size(MSA,2)))
        setannotfile!(annot, "ColMap", join(vcat(1:size(MSA,2)), ','))
        for i in 1:length(IDS)
            setannotsequence!(annot, IDS[i], "SeqMap", MAP[i])
        end
    else
        MSA = NamedArray(convert(Matrix{Residue}, SEQS))
        setnames!(MSA, IDS, 1)
    end
    msa = AnnotatedMultipleSequenceAlignment(MSA, annot)
    if deletefullgaps
        deletefullgapcolumns!(msa)
    end
    msa
end

function Base.parse(io::Union{IO, AbstractString},
                   format::Type{Stockholm},
                   output::Type{MultipleSequenceAlignment};
                   deletefullgaps::Bool=true)
    IDS, SEQS = _pre_readstockholm_sequences(io)
    msa = MultipleSequenceAlignment(convert(Matrix{Residue}, SEQS))
    setnames!(namedmatrix(msa), IDS, 1)
    if deletefullgaps
        deletefullgapcolumns!(msa)
    end
    msa
end

function Base.parse(io::Union{IO,AbstractString},
                   format::Type{Stockholm},
                   output::Type{Matrix{Residue}};
                   deletefullgaps::Bool=true)
    IDS, SEQS = _pre_readstockholm_sequences(io)
    _strings_to_msa(SEQS, deletefullgaps, checkalphabet)
end

function Base.parse(io, format::Type{Stockholm};
                    generatemapping::Bool=false,
                    useidcoordinates::Bool=false,
                    deletefullgaps::Bool=true,
                    keepinserts::Bool=false)
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
        if seq_id in keys(seq_dict)
            push!(seq_dict[seq_id], string(key[2], '\t', value))
        else
            seq_dict[seq_id] = [ string(key[2], '\t', value) ]
        end
    end
    sizehint!(seq_dict, length(seq_dict))
end

function Base.print(io::IO, msa::AnnotatedMultipleSequenceAlignment, format::Type{Stockholm})
    _printfileannotations(io, msa.annotations)
    _printsequencesannotations(io, msa.annotations)
    res_annotations = _to_sequence_dict(msa.annotations.residues)
    for i in 1:nsequences(msa)
        id = msa.id[i]
        seq = stringsequence(msa, i)
        println(io, string(id, "\t\t", seq))
        if id in keys(res_annotations)
            for line in res_annotations[id]
                println(io, string("#=GR ", id, '\t', line))
            end
        end
    end
    _printcolumnsannotations(io, msa.annotations)
    println(io, "//")
end

function Base.print(io::IO, msa::MultipleSequenceAlignment, format::Type{Stockholm})
    for i in 1:nsequences(msa)
        id = msa.id[i]
        seq = stringsequence(msa, i)
        println(io, string(id, "\t\t", seq))
    end
    println(io, "//")
end

Base.print(msa::AnnotatedMultipleSequenceAlignment) = print(STDOUT, msa, Stockholm)
