immutable Stockholm <: Format end

function _fill_with_line!(line, IDS, SEQS, GF, GS, GC, GR)
    if line[1:4] == "#=GF"
        words = get_n_words(line,3)
        id = words[2]
        if id in keys(GF)
            GF[ id ] = GF[ id ] * "\n" * words[3]
        else
            GF[ id ] = words[3]
        end
    elseif line[1:4] == "#=GS"
        words = get_n_words(line,4)
        idtuple = (words[2],words[3])
        if idtuple in keys(GS)
            GS[ idtuple ] = GS[ idtuple ] * "\n" * words[4]
        else
            GS[ idtuple ] = words[4]
        end
    elseif line[1:4] == "#=GC"
        words = get_n_words(line,3)
        GC[words[2]] = words[3]
    elseif line[1:4] == "#=GR"
        words = get_n_words(line,4)
        GR[(words[2],words[3])] = words[4]
    elseif line[1:1] != "#"
        words = get_n_words(line,2)
        push!(IDS, words[1])
        push!(SEQS,words[2])
    end
end

function _pre_readstockholm(io::Union{IO, AbstractString})
    IDS  = ASCIIString[]
    SEQS = ASCIIString[]
    GF = OrderedDict{ASCIIString,ByteString}()
    GC = Dict{ASCIIString,ASCIIString}()
    GS = Dict{Tuple{ASCIIString,ASCIIString},ASCIIString}()
    GR = Dict{Tuple{ASCIIString,ASCIIString},ASCIIString}()

    for line in eachline(io)
        if length(line) >= 4
            _fill_with_line!(line, IDS, SEQS, GF, GS, GC, GR)
        end
    end

    #GF = sizehint!(GF, length(GF))
    GF = sizehint!(GF, length(GF))
    GC = sizehint!(GC, length(GC))
    GS = sizehint!(GS, length(GS))
    GR = sizehint!(GR, length(GR))
    (IDS, SEQS, GF, GS, GC, GR)
end

function parse(io::Union{IO, AbstractString}, format::Type{Stockholm},
               output::Type{AnnotatedMultipleSequenceAlignment}; generatemapping::Bool=false,
               useidcoordinates::Bool=false, deletefullgaps::Bool=true, checkalphabet::Bool=false)
    IDS, SEQS, GF, GS, GC, GR = _pre_readstockholm(io)
    annot = Annotations(GF, GS, GC, GR)
    if generatemapping
        MSA, MAP = useidcoordinates && hascoordinates(IDS[1]) ? _to_msa_mapping(SEQS, IDS) : _to_msa_mapping(SEQS)
        setannotfile!(annot, "NCol", string(size(MSA,2)))
        setannotfile!(annot, "ColMap", join(vcat(1:size(MSA,2)), ','))
        for i in 1:length(IDS)
            setannotsequence!(annot, IDS[i], "SeqMap", MAP[i])
        end
    else
        MSA = convert(Matrix{Residue}, SEQS)
    end
    msa = AnnotatedMultipleSequenceAlignment(IndexedArray(IDS), MSA, annot)
    if checkalphabet
        deletenotalphabetsequences!(msa, SEQS)
    end
    if deletefullgaps
        deletefullgapcolumns!(msa)
    end
    msa
end

function parse(io::Union{IO, AbstractString}, format::Type{Stockholm}, output::Type{MultipleSequenceAlignment}; deletefullgaps::Bool=true, checkalphabet::Bool=false)
    # Could be faster with a special _pre_readstockholm
    IDS, SEQS, GF, GS, GC, GR = _pre_readstockholm(io)
    msa = MultipleSequenceAlignment(IndexedArray(IDS), convert(Matrix{Residue}, SEQS))
    if checkalphabet
        deletenotalphabetsequences!(msa, SEQS)
    end
    if deletefullgaps
        deletefullgapcolumns!(msa)
    end
    msa
end

function parse(io::Union{IO,AbstractString}, format::Type{Stockholm}, output::Type{Matrix{Residue}}; deletefullgaps::Bool=true, checkalphabet::Bool=false)
    # Could be faster with a special _pre_readstockholm
    IDS, SEQS, GF, GS, GC, GR = _pre_readstockholm(io)
    _strings_to_msa(SEQS, deletefullgaps, checkalphabet)
end

parse(io, format::Type{Stockholm};  generatemapping::Bool=false,
      useidcoordinates::Bool=false, deletefullgaps::Bool=true, checkalphabet::Bool=false) = parse(io, Stockholm, AnnotatedMultipleSequenceAlignment,
                                                                                                  generatemapping=generatemapping,
                                                                                                  useidcoordinates=useidcoordinates,
                                                                                                  deletefullgaps=deletefullgaps,
                                                                                                  checkalphabet=checkalphabet)

# Print Pfam
# ==========

function _to_sequence_dict(annotation::Dict{Tuple{ASCIIString,ASCIIString},ASCIIString})
    seq_dict = Dict{ASCIIString,Vector{ASCIIString}}()
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

function print(io::IO, msa::AnnotatedMultipleSequenceAlignment, format::Type{Stockholm})
    _printfileannotations(io, msa.annotations)
    _printsequencesannotations(io, msa.annotations)
    res_annotations = _to_sequence_dict(msa.annotations.residues)
    for i in 1:nsequences(msa)
        id = msa.id[i]
        seq = asciisequence(msa, i)
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

function print(io::IO, msa::MultipleSequenceAlignment, format::Type{Stockholm})
    for i in 1:nsequences(msa)
        id = msa.id[i]
        seq = asciisequence(msa, i)
        println(io, string(id, "\t\t", seq))
    end
    println(io, "//")
end

print(msa::AnnotatedMultipleSequenceAlignment) = print(STDOUT, msa, Stockholm)
