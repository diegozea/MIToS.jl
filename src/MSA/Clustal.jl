struct Clustal <: MSAFormat end

# Parser based on the format description at
# https://meme-suite.org/meme/doc/clustalw-format.html
# Each block ends with a line showing the degree of conservation for
# the columns of the alignment.

# Match a header token, not a sequence name such as CLUSTAL_seq.
"""
Return the regular expression that identifies an alignment header.
Used by `hasnextmsa` to check for another alignment and by parsers to recognize
where the next alignment begins, so they can stop reading the current one.
"""
_msa_header(::Type{Clustal}) = r"^CLUSTALW?(?:\s|$)"

"""
Read one Clustal alignment. Put the next header back in buffered inputs;
otherwise consume it as with sequential line reading.
Return sequence names, sequences and annotations.
"""
function _load_clustal_sequences(io::IO)
    seqs = OrderedDict{String,String}()
    conservation = IOBuffer()
    seq_re = r"^(\S+)\s+([A-Za-z.-]+)(?:\s+\d+)?"  # sequence line with optional count
    startidx = 0
    endidx = 0
    seen_header = false
    in_sequence_block = false # true when reading a sequence block
    buffered = io isa TranscodingStream
    while !eof(io)
        line = _read_msa_line(io)
        chomped = chomp(line)
        stripped = strip(chomped)
        # blank line ends the current sequence block
        isempty(stripped) && (in_sequence_block = false; continue)
        if occursin(_msa_header(Clustal), stripped)
            if seen_header || !isempty(seqs)
                buffered && TranscodingStreams.unread(io, codeunits(line))
                break
            end
            seen_header = true
            continue
        end
        startswith(chomped, '#') && continue
        if (m = match(seq_re, chomped)) !== nothing  # sequence line
            id = m.captures[1]
            seq = m.captures[2]
            if !in_sequence_block
                # store positions of the sequence slice to read conservation line
                startidx = m.offsets[2]
                endidx = startidx + length(seq) - 1
            end
            if haskey(seqs, id)
                seqs[id] = seqs[id] * seq
            else
                seqs[id] = seq
            end
            in_sequence_block = true  # we are inside a sequence block now
            continue
        end
        # Preserve conservation spaces, but exclude line endings.
        if in_sequence_block && isascii(chomped)
            # conservation line found
            stop = min(endidx, lastindex(chomped))
            if stop >= startidx
                # remove leading/trailing padding spaces from the conservation
                # line before storing it using the previously stored indices
                # of the aligned columns in this block.
                consblock = chomped[startidx:stop]
                write(conservation, consblock)
            end
            in_sequence_block = false
        end
    end
    IDS = collect(keys(seqs))
    SEQS = collect(values(seqs))
    CONS = String(take!(conservation))
    annot = Annotations()
    _disambiguate_seqnames!(IDS, annot)
    isempty(CONS) || setannotcolumn!(annot, "cons", CONS)
    return IDS, SEQS, annot
end

function _load_sequences(
    io::Union{IO,AbstractString},
    format::Type{Clustal};
    create_annotations::Bool = false,
)
    _load_clustal_sequences(io isa AbstractString ? IOBuffer(io) : io)
end

function Utils.print_file(
    io::IO,
    msa::AbstractMatrix{Residue},
    format::Type{Clustal};
    showcounts::Bool = false,
)
    seqnames = sequencenames(msa)
    namew = maximum(length.(seqnames)) + 2
    println(io, "CLUSTAL\n")
    ncol = ncolumns(msa)
    block = 60
    cons = nothing
    if isa(msa, AnnotatedAlignedObject)
        cons = getannotcolumn(msa, "cons", "")
    end
    # Clustal prints blocks of 60 columns
    for start = 1:block:ncol  # start of a new block
        stop = min(start + block - 1, ncol)
        for i = 1:nsequences(msa)
            seq = stringsequence(getsequence(msa, i)[:, start:stop])
            line = rpad(seqnames[i], namew) * seq
            if showcounts
                # append cumulative residue count as Clustal does
                pre = stringsequence(getsequence(msa, i)[:, 1:stop])
                rescount = Base.count(ch -> ch != '-' && ch != '.', pre)
                line *= " " * string(rescount)
            end
            println(io, line)
        end
        if cons !== nothing
            println(io, rpad("", namew), cons[start:stop])
        end
        println(io)  # blank line between blocks
    end
end
