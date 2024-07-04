struct Stockholm <: MSAFormat end

@inline function _fill_with_sequence_line!(IDS, SEQS, line)
    if !startswith(line, '#') && !startswith(line, "//")
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
    if startswith(line, "#=GF")
        words = get_n_words(line, 3)
        id = words[2]
        if id in keys(GF)
            GF[id] = GF[id] * "\n" * words[3]
        else
            GF[id] = words[3]
        end
    elseif startswith(line, "#=GS")
        words = get_n_words(line, 4)
        idtuple = (words[2], words[3])
        if idtuple in keys(GS)
            GS[idtuple] = GS[idtuple] * "\n" * words[4]
        else
            GS[idtuple] = words[4]
        end
    elseif startswith(line, "#=GC")
        words = get_n_words(line, 3)
        GC[words[2]] = words[3]
    elseif startswith(line, "#=GR")
        words = get_n_words(line, 4)
        GR[(words[2], words[3])] = words[4]
    else
        _fill_with_sequence_line!(IDS, SEQS, line)
    end
end

function _pre_readstockholm(io::Union{IO,AbstractString})
    IDS = OrderedSet{String}()
    SEQS = String[]
    GF = OrderedDict{String,String}()
    GC = Dict{String,String}()
    GS = Dict{Tuple{String,String},String}()
    GR = Dict{Tuple{String,String},String}()

    @inbounds for line::String in lineiterator(io)
        isempty(line) && continue
        startswith(line, "//") && break
        _fill_with_line!(IDS, SEQS, GF, GS, GC, GR, line)
    end

    GF = sizehint!(GF, length(GF))
    GC = sizehint!(GC, length(GC))
    GS = sizehint!(GS, length(GS))
    GR = sizehint!(GR, length(GR))
    (IDS, SEQS, GF, GS, GC, GR)
end

function _pre_readstockholm_sequences(io::Union{IO,AbstractString})
    IDS = OrderedSet{String}()
    SEQS = String[]
    @inbounds for line::String in lineiterator(io)
        isempty(line) && continue
        startswith(line, "//") && break
        _fill_with_sequence_line!(IDS, SEQS, line)
    end
    (IDS, SEQS)
end

function _load_sequences(
    io::Union{IO,AbstractString},
    format::Type{Stockholm};
    create_annotations::Bool = false,
)
    if create_annotations
        IDS, SEQS, GF, GS, GC, GR = _pre_readstockholm(io)
        annot = Annotations(GF, GS, GC, GR)
    else
        IDS, SEQS = _pre_readstockholm_sequences(io)
        annot = Annotations()
    end
    return collect(IDS), SEQS, annot
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
            seq_dict[seq_id] = [string(seq_id, '\t', key[2], '\t', value)]
        end
    end
    sizehint!(seq_dict, length(seq_dict))
end

function Utils.print_file(io::IO, msa::AbstractMatrix{Residue}, format::Type{Stockholm})
    has_annotations = isa(msa, AnnotatedAlignedObject) && !isempty(msa.annotations)
    if has_annotations
        _printfileannotations(io, msa.annotations)
        _printsequencesannotations(io, msa.annotations)
        res_annotations = _to_sequence_dict(msa.annotations.residues)
    end
    seqnames = sequencenames(msa)
    aligned = _get_aligned_columns(msa)
    for i = 1:nsequences(msa)
        id = seqnames[i]
        seq = stringsequence(msa, i)
        formatted_seq = _format_inserts(seq, aligned)
        println(io, id, "\t\t\t", formatted_seq)
        if has_annotations && haskey(res_annotations, id)
            for line in res_annotations[id]
                println(io, "#=GR ", line)
            end
        end
    end
    has_annotations && _printcolumnsannotations(io, msa.annotations)
    println(io, "//")
end
