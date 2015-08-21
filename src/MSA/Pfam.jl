using MIToS.Utils

import Base: parse, print

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

function _pre_readstockholm(io::Union(IO, AbstractString))
  IDS  = ASCIIString[]
  SEQS = ASCIIString[]
  GF = Dict{ASCIIString,ASCIIString}()
  GC = Dict{ASCIIString,ASCIIString}()
  GS = Dict{Tuple{ASCIIString,ASCIIString},ASCIIString}()
  GR = Dict{Tuple{ASCIIString,ASCIIString},ASCIIString}()

  for line in eachline(io)
     if length(line) >= 4
       _fill_with_line!(line, IDS, SEQS, GF, GS, GC, GR)
     end
  end

  GF = sizehint!(GF, length(GF))
  GC = sizehint!(GC, length(GC))
  GS = sizehint!(GS, length(GS))
  GR = sizehint!(GR, length(GR))
  (IDS, SEQS, GF, GS, GC, GR)
end

function parse(io::Union(IO, AbstractString), format::Type{Stockholm}, output::Type{AnnotatedMultipleSequenceAlignment}; useidcoordinates::Bool=true, deletefullgaps::Bool=true)
  IDS, SEQS, GF, GS, GC, GR = _pre_readstockholm(io)
  MSA, MAP = useidcoordinates && hascoordinates(IDS[1]) ? _to_msa_mapping(SEQS, IDS) : _to_msa_mapping(SEQS)
  COLS = vcat(1:size(MSA,2))
  msa = AnnotatedMultipleSequenceAlignment(IndexedVector(IDS), MSA, MAP, IndexedVector(COLS), Annotations(GF, GS, GC, GR))
  if deletefullgaps
    deletefullgaps!(msa)
  end
  return(msa)
end

function parse(io::Union(IO, AbstractString), format::Type{Stockholm}, output::Type{MultipleSequenceAlignment}; deletefullgaps::Bool=true)
  # Could be faster with a special _pre_readstockholm
  IDS, SEQS, GF, GS, GC, GR = _pre_readstockholm(io)
  msa = MultipleSequenceAlignment(IndexedVector(IDS), convert(Matrix{Residue}, SEQS))
  if deletefullgaps
    deletefullgaps!(msa)
  end
  return(msa)
end

parse(io, format::Type{Stockholm}; useidcoordinates::Bool=true, deletefullgaps::Bool=true) = parse(io, Stockholm, AnnotatedMultipleSequenceAlignment, useidcoordinates=useidcoordinates, deletefullgaps=deletefullgaps)

# Print Pfam
# ==========

function _printfileannotations(io::IO, msa::AnnotatedMultipleSequenceAlignment)
	if !isempty(msa.annotations.file)
		for (key, value) in msa.annotations.file
			println(io, string("#=GF ", key, " ", value))
		end
	end
end

function _printcolumnsannotations(io::IO, msa::AnnotatedMultipleSequenceAlignment)
	if !isempty(msa.annotations.columns)
		for (key, value) in msa.annotations.columns
			println(io, string("#=GC ", key, "\t\t", value))
		end
	end
end

function _printsequencesannotations(io::IO, msa::AnnotatedMultipleSequenceAlignment)
	if !isempty(msa.annotations.sequences)
		for (key, value) in msa.annotations.sequences
			println(io, string("#=GS ", key[1], " ", key[2], " ", value))
		end
	end
end


function _to_sequence_dict(annotation::Dict{Tuple{ASCIIString,ASCIIString},ASCIIString})
	seq_dict = Dict{ASCIIString,Vector{ASCIIString}}()
	for (key, value) in annotation
		seq_id = key[1]
		if seq_id in keys(seq_dict)
			push!(seq_dict[seq_id], string(key[2], "\t", value))
		else
			seq_dict[seq_id] = [ string(key[2], "\t", value) ]
		end
	end
	sizehint!(seq_dict, length(seq_dict))
end

function print(io::IO, msa::AnnotatedMultipleSequenceAlignment, format::Type{Stockholm})
	_printfileannotations(io, msa)
	_printsequencesannotations(io, msa)
	res_annotations = _to_sequence_dict(msa.annotations.residues)
	for i in 1:nsequences(msa)
		id = selectvalue(msa.id, i)
		seq = asciisequence(msa, i)
		println(io, string(id, "\t\t", seq))
		if id in keys(res_annotations)
			for line in res_annotations[id]
				println(io, string("#=GR ", id, " ", line))
			end
		end
	end
	_printcolumnsannotations(io, msa)
	println(io, "//")
end

function print(io::IO, msa::MultipleSequenceAlignment, format::Type{Stockholm})
	for i in 1:nsequences(msa)
		id = selectvalue(msa.id, i)
		seq = asciisequence(msa, i)
		println(io, string(id, "\t\t", seq))
	end
	println(io, "//")
end

print(msa::AnnotatedMultipleSequenceAlignment) = print(STDOUT, msa, Stockholm)

# Download Pfam
# =============
"""
Download a gzipped stockholm full alignment for the `pfamcode`.
The extension of the downloaded file is `.stockholm.gz` by default. The `filename` can be changed, but the `.gz` at the end is mandatory.
"""
function downloadpfam(pfamcode::ASCIIString; filename::ASCIIString="$pfamcode.stockholm.gz")
  if ismatch(r"^PF\d{5}$"i, pfamcode) # length(pfamcode)== 7 && ( pfamcode[1:2] == "PF" || pfamcode[1:2] == "pf" )
    #namegz = string(filename, ".gz")
    download(string("http://pfam.xfam.org/family/PF", pfamcode[3:end], "/alignment/full/gzipped"), filename) # namegz)
    #run(`gzip -d $namegz`)
  else
    throw( ErrorException( string(pfamcode, " is not a correct Pfam code") ) )
  end
end
