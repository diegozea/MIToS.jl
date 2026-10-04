struct Clustal <: MSAFormat end

# Parser based on the format description at
# https://meme-suite.org/meme/doc/clustalw-format.html
# Each block ends with a line showing the degree of conservation for
# the columns of the alignment.

# Match a header token, not a sequence name such as CLUSTAL_seq.
_msa_header(::Type{Clustal}) = r"^CLUSTALW?(?:\s|$)"

"""
Read sequence data and conservation annotations from an iterable of Clustal lines.
Return sequence names, sequences, annotations and whether another header was consumed.
Set `header_read` when the caller has already consumed the current alignment's header.
"""
function _load_clustal_sequences(
    lines;
    header_read::Bool = false,
    fail_on_duplicate_seqnames::Bool = false,
)
    seqs = OrderedDict{String,String}()
    block_ids = fail_on_duplicate_seqnames ? Set{String}() : nothing
    conservation = IOBuffer()
    seq_re = r"^(\S+)\s+([A-Za-z.-]+)(?:\s+\d+)?"  # sequence line with optional count
    startidx = 0
    endidx = 0
    seen_header = header_read
    has_next = false
    in_sequence_block = false # true when reading a sequence block
    for line in lines
        chomped = chomp(line)
        # blank line ends the current sequence block
        isempty(strip(chomped)) && (in_sequence_block = false; continue)
        if occursin(_msa_header(Clustal), chomped)
            # A new header starts another alignment, not another sequence block.
            if seen_header || !isempty(seqs)
                has_next = true
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
                fail_on_duplicate_seqnames && empty!(block_ids)
                # store positions of the sequence slice to read conservation line
                startidx = m.offsets[2]
                endidx = startidx + length(seq) - 1
            end
            fail_on_duplicate_seqnames && _check_unique_seqname!(block_ids, id)
            if haskey(seqs, id)
                seqs[id] = seqs[id] * seq
            else
                seqs[id] = seq
            end
            in_sequence_block = true  # we are inside a sequence block now
            continue
        end
        # using line instead of chomped to preserve whitespaces in the conservation line
        if in_sequence_block && isascii(line) && match(seq_re, line) === nothing
            # conservation line found
            stop = min(endidx, lastindex(line))
            if stop >= startidx
                # remove leading/trailing padding spaces from the conservation
                # line before storing it using the previously stored indices
                # of the aligned columns in this block.
                consblock = line[startidx:stop]
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
    return IDS, SEQS, annot, has_next
end

function _load_sequences(
    io::Union{IO,AbstractString},
    format::Type{Clustal};
    create_annotations::Bool = false,
    fail_on_duplicate_seqnames::Bool = false,
)
    _load_clustal_sequences(
        lineiterator(io);
        fail_on_duplicate_seqnames = fail_on_duplicate_seqnames,
    )
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
