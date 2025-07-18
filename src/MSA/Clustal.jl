struct Clustal <: MSAFormat end

# Parser based on the format description at
# https://meme-suite.org/meme/doc/clustalw-format.html
# Each block ends with a line showing the degree of conservation for
# the columns of the alignment.

function _pre_readclustal(io::Union{IO,AbstractString})
    seqs = OrderedDict{String,String}()
    conservation = IOBuffer()
    seq_re = r"^(\S+)\s+([A-Za-z.-]+)(?:\s+\d+)?"  # sequence line with optional count
    startidx = 0
    endidx = 0
    in_sequence_block = false # true when reading a sequence block
    for line in lineiterator(io)
        chomped = chomp(line)
        # blank line ends the current sequence block
        isempty(strip(chomped)) && (in_sequence_block = false; continue)
        startswith(chomped, "CLUSTAL") && continue
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
        if in_sequence_block && isascii(chomped) && match(seq_re, chomped) === nothing
            # conservation line found
            stop = min(endidx, lastindex(chomped))
            if stop >= startidx
                # remove leading/trailing padding spaces from the conservation
                # line before storing it. The remaining characters correspond
                # to the alignment columns in this block.
                consblock = strip(chomped[startidx:stop])
                write(conservation, consblock)
            end
            in_sequence_block = false
        end
    end
    IDS = collect(keys(seqs))
    SEQS = collect(values(seqs))
    CONS = String(take!(conservation))
    (IDS, SEQS, isempty(CONS) ? nothing : CONS)
end

function _load_sequences(
    io::Union{IO,AbstractString},
    format::Type{Clustal};
    create_annotations::Bool = false,
)
    IDS, SEQS, CONS = _pre_readclustal(io)
    annot = Annotations()
    _disambiguate_seqnames!(IDS, annot)
    CONS !== nothing && setannotcolumn!(annot, "cons", CONS)
    return IDS, SEQS, annot
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
