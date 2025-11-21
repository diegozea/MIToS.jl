struct Stockholm <: MSAFormat end

# NOTE: Sequence‑Name Disambiguation
# We do not support sequence‑name disambiguation via the `OnlineSequenceNameDisambiguator`
# in Stockholm format, because duplicate sequence names are not permitted.
# Long sequences are split across multiple blocks, and the name is used to
# reassemble the fragments and to key per‑sequence annotations. An example of
# file with multiple blocks is at test/data/clustalo-I20240512-trunc.aln-stockholm

@inline function _fill_with_sequence_line!(ids, seqs, id_to_index, line)
    if !startswith(line, '#') && !startswith(line, "//")
        words = get_n_words(line, 2)
        @inbounds id = words[1]
        idx = get(id_to_index, id, 0)
        if idx == 0
            push!(ids, id)
            push!(seqs, words[2])
            id_to_index[id] = length(seqs)
            return
        end
        seq = seqs[idx]
        if seq isa String
            buffer = IOBuffer(; read = true, write = true)
            write(buffer, seq)
            write(buffer, words[2])
            seqs[idx] = buffer
        else
            write(seq, words[2])
        end
    end
end

@inline function _append_annotation_fragment!(dict::AbstractDict, key, fragment)
    if haskey(dict, key)
        dict[key] = dict[key] * fragment
    else
        dict[key] = fragment
    end
end

function _fill_with_line!(IDS, SEQS, id_to_index, GF, GS, GC, GR, line)
    if startswith(line, "#=GF")
        words = get_n_words(line, 3)
        id = words[2]
        if haskey(GF, id)
            GF[id] = GF[id] * "\n" * words[3]
        else
            GF[id] = words[3]
        end
    elseif startswith(line, "#=GS")
        words = get_n_words(line, 4)
        idtuple = (words[2], words[3])
        if haskey(GS, idtuple)
            GS[idtuple] = GS[idtuple] * "\n" * words[4]
        else
            GS[idtuple] = words[4]
        end
    elseif startswith(line, "#=GC")
        words = get_n_words(line, 3)
        _append_annotation_fragment!(GC, words[2], words[3])
    elseif startswith(line, "#=GR")
        words = get_n_words(line, 4)
        key = (words[2], words[3])
        fragment = words[4]
        _append_annotation_fragment!(GR, key, fragment)
    else
        _fill_with_sequence_line!(IDS, SEQS, id_to_index, line)
    end
end

function _pre_readstockholm(io::Union{IO,AbstractString})
    IDS = String[]
    seqs = Union{String,IOBuffer}[]
    id_to_index = Dict{String,Int}()
    GF = OrderedDict{String,String}()
    GC = Dict{String,String}()
    GS = Dict{Tuple{String,String},String}()
    GR = Dict{Tuple{String,String},String}()

    @inbounds for line in lineiterator(io)
        isempty(line) && continue
        startswith(line, "//") && break
        _fill_with_line!(IDS, seqs, id_to_index, GF, GS, GC, GR, line)
    end

    GF = sizehint!(GF, length(GF))
    GC = sizehint!(GC, length(GC))
    GS = sizehint!(GS, length(GS))
    GR = sizehint!(GR, length(GR))
    (IDS, [seq isa String ? seq : String(take!(seq)) for seq in seqs], GF, GS, GC, GR)
end

function _pre_readstockholm_sequences(io::Union{IO,AbstractString})
    IDS = String[]
    seqs = Union{String,IOBuffer}[]
    id_to_index = Dict{String,Int}()
    @inbounds for line in lineiterator(io)
        isempty(line) && continue
        startswith(line, "//") && break
        _fill_with_sequence_line!(IDS, seqs, id_to_index, line)
    end
    (IDS, [seq isa String ? seq : String(take!(seq)) for seq in seqs])
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
    return IDS, SEQS, annot
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
