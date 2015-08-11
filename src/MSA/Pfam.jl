using MIToS.Utils

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

function _pre_readstockholm(filename::ASCIIString)
  IDS  = ASCIIString[]
  SEQS = ASCIIString[]
  GF = Dict{ASCIIString,ASCIIString}()
  GC = Dict{ASCIIString,ASCIIString}()
  GS = Dict{Tuple{ASCIIString,ASCIIString},ASCIIString}()
  GR = Dict{Tuple{ASCIIString,ASCIIString},ASCIIString}()

  open(filename, "r") do stockholm_file
    for line in eachline(stockholm_file)
      if length(line) >= 4
        _fill_with_line!(line, IDS, SEQS, GF, GS, GC, GR)
     end
    end
  end

  GF = sizehint!(GF, length(GF))
  GC = sizehint!(GC, length(GC))
  GS = sizehint!(GS, length(GS))
  GR = sizehint!(GR, length(GR))
  (IDS, SEQS, GF, GS, GC, GR)
end

function _to_msa_mapping(sequences::Array{ASCIIString,1})
  nseq = size(sequences,1)
  nres = length(sequences[1])
  aln = Array(Residue,nres,nseq)
  mapp = zeros(Int,nseq,nres) # This needs to be zeros
  gaps = UInt8['.', '-']
  for i in 1:nseq
    init = 1
    seq = sequences[i].data
    if length(seq) == nres
      @inbounds for j in 1:nres
        res = seq[j]
        aln[j,i] = Residue( res )
        if !( res in gaps )
          mapp[i, j] = init
          init += 1
        end
      end
    else
      throw( ErrorException( "There is and aligned sequence with different number of columns [ $(length(seq)) != $(nres) ]:\n$(ascii(seq))" ) )
    end
  end
  (aln', mapp)
end

function _to_msa_mapping(sequences::Array{ASCIIString,1}, ids::Array{ASCIIString,1})
  nseq = size(sequences,1)
  nres = length(sequences[1])
  aln = Array(Residue,nres,nseq)
  mapp = zeros(Int,nseq,nres) # This needs to be zeros
  gaps = UInt8['.', '-']
  sep = r"/|-"
  for i in 1:nseq
    fields = split(ids[i],sep)
    init = length(fields) == 3 ? parse(Int, fields[2]) : 1
    seq = sequences[i].data
    if length(seq) == nres
      @inbounds for j in 1:nres
        res = seq[j]
        aln[j,i] = Residue( res )
        if !( res in gaps )
          mapp[i, j] = init
          init += 1
        end
      end
    else
      throw( ErrorException( "There is and aligned sequence with different number of columns [ $(length(seq)) != $(nres) ]:\n$(ascii(seq))" ) )
    end
    if (init - 1) != parse(Int, fields[3])
      throw( ErrorException( "Different lengths: $(fields[3]) != $(init) for sequence $(ids[i])" ) )
    end
  end
  (aln', mapp)
end

function readpfam(filename::ASCIIString; useidcoordinates::Bool=true)
  IDS, SEQS, GF, GS, GC, GR = _pre_readstockholm(filename)
  MSA, MAP = useidcoordinates && hascoordinates(IDS[1]) ? _to_msa_mapping(SEQS, IDS) : _to_msa_mapping(SEQS)
  COLS = vcat(1:size(MSA,2))
  msa = MultipleSequenceAlignment(IndexedVector(IDS), MSA, MAP, IndexedVector(COLS), Annotations(GF, GS, GC, GR))
  filtercolumns!(msa, columngappercentage(msa) .< 1.0)
end

# Print Pfam
# ==========

function _printfileannotations(io::IO, msa::MultipleSequenceAlignment)
	if !isempty(msa.annotations.file)
		for (key, value) in msa.annotations.file
			println(io, string("#=GF ", key, " ", value))
		end
	end
end

function _printcolumnsannotations(io::IO, msa::MultipleSequenceAlignment)
	if !isempty(msa.annotations.columns)
		for (key, value) in msa.annotations.columns
			println(io, string("#=GC ", key, "\t\t", value))
		end
	end
end

function _printsequencesannotations(io::IO, msa::MultipleSequenceAlignment)
	if !isempty(msa.annotations.sequences)
		for (key, value) in msa.annotations.sequences
			println(io, string("#=GC ", key[1], " ", key[2], " ", value))
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

function printpfam(io::IO, msa::MultipleSequenceAlignment)
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

printpfam(msa::MultipleSequenceAlignment) = printpfam(STDOUT, msa)

# Write Pfam
# ==========

function writepfam(filename::ASCIIString, msa::MultipleSequenceAlignment)
	open(filename, "w") do fh
		printpfam(fh, msa)
	end
end

# Download Pfam
# =============

function downloadpfam(pfamcode::ASCIIString; filename::ASCIIString="$pfamcode.aln")
  if length(pfamcode)== 7 && ( pfamcode[1:2] == "PF" || pfamcode[1:2] == "pf" )
    namegz = string(filename, ".gz")
    download(string("http://pfam.xfam.org/family/PF", pfamcode[3:end], "/alignment/full/gzipped"), namegz)
    run(`gzip -d $namegz`)
  else
    throw( ErrorException( string(pfamcode, " is not a correct Pfam code") ) )
  end
end
