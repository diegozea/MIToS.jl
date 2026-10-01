struct Clustal <: MSAFormat end

# Parser based on the format description at
# https://meme-suite.org/meme/doc/clustalw-format.html
# Each block ends with a line showing the degree of conservation for
# the columns of the alignment.

# Match a header token, not a sequence name such as CLUSTAL_seq.
_is_clustal_header(line::AbstractString) = occursin(r"^CLUSTALW?(?:\s|$)", line)

_pre_readclustal(io::Union{IO,AbstractString}) = _pre_readclustal(lineiterator(io))

function _pre_readclustal(lines)
    seqs = OrderedDict{String,String}()
    conservation = IOBuffer()
    seq_re = r"^(\S+)\s+([A-Za-z.-]+)(?:\s+\d+)?"  # sequence line with optional count
    startidx = 0
    endidx = 0
    in_sequence_block = false # true when reading a sequence block
    for line in lines
        chomped = chomp(line)
        # blank line ends the current sequence block
        isempty(strip(chomped)) && (in_sequence_block = false; continue)
        _is_clustal_header(chomped) && continue
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
    (IDS, SEQS, isempty(CONS) ? nothing : CONS)
end

"""
Read sequence data and conservation annotations from an iterable of Clustal lines.
"""
function _load_clustal_sequences(lines)
    IDS, SEQS, CONS = _pre_readclustal(lines)
    annot = Annotations()
    _disambiguate_seqnames!(IDS, annot)
    CONS !== nothing && setannotcolumn!(annot, "cons", CONS)
    return IDS, SEQS, annot
end

function _load_sequences(
    io::Union{IO,AbstractString},
    format::Type{Clustal};
    create_annotations::Bool = false,
)
    _load_clustal_sequences(lineiterator(io))
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
