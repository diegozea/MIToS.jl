function __words(line::ASCIIString, n::Int)
  N = length(line.data)
  coords = Array(ASCIIString, n)
  spaces = Uint8[' ', '\n', '\t']
  j = 0
  start = 0
  if N > 1
    block = !(line.data[1] in spaces)
    if block
      start = 1
    end
    for i in 2:N
      if line.data[i] in spaces
        if block
          j += 1
          if j == n
            last = line.data[N] == uint8('\n') ? N-1 : N
            @inbounds coords[j] = ascii(getindex(line.data,start:last))
            break
          else
            @inbounds coords[j] = ascii(getindex(line.data,start:i-1))
          end
        end
        block = false
      else
        if ! block
          block = true
          start = i
        end
      end
    end
  elseif N==1
    if line.data[1] in spaces
      j=0
    else
      j=1
      coords[1] = line
    end
  else
    j=0
  end
  if j != n
    resize!(coords, j)
  end
  coords
end

function __fill_with_line!(line, IDS, SEQS, GF, GS, GC, GR)
  if line[1:4] == "#=GF"
    words = __words(line,3)
    id = words[2]
    if id in keys(GF)
      GF[ id ] = GF[ id ] * "\n" * words[3]
    else
      GF[ id ] = words[3]
    end
  elseif line[1:4] == "#=GS"
    words = __words(line,4)
    idtuple = (words[2],words[3])
    if idtuple in keys(GS)
      GS[ idtuple ] = GS[ idtuple ] * "\n" * words[4]
    else
      GS[ idtuple ] = words[4]
    end
  elseif line[1:4] == "#=GC"
    words = __words(line,3)
    GC[words[2]] = words[3]
  elseif line[1:4] == "#=GR"
    words = __words(line,4)
    GR[(words[2],words[3])] = words[4]
  elseif line[1:1] != "#"
    words = __words(line,2)
    push!(IDS, words[1])
    push!(SEQS,words[2])
  end
end

function __pre_readstockholm(filename::ASCIIString)
  IDS  = ASCIIString[]
  SEQS = ASCIIString[]
  GF = Dict{ASCIIString,ASCIIString}()
  GC = Dict{ASCIIString,ASCIIString}()
  GS = Dict{(ASCIIString,ASCIIString),ASCIIString}()
  GR = Dict{(ASCIIString,ASCIIString),ASCIIString}()

  open(filename, "r") do stockholm_file
    for line in eachline(stockholm_file)
      if length(line) >= 4
        __fill_with_line!(line, IDS, SEQS, GF, GS, GC, GR)
     end
    end
  end

  GF = sizehint(GF, length(GF))
  GC = sizehint(GC, length(GC))
  GS = sizehint(GS, length(GS))
  GR = sizehint(GR, length(GR))
  (IDS, SEQS, GF, GS, GC, GR)
end

function __to_msa(sequences::Array{ASCIIString,1})
  nseq = length(sequences)
  nres = length(sequences[1])
  aln = Array(Residue,nres,nseq)
  for col in 1:nseq
    seq = sequences[col].data
    if length(seq) == nres
      aln[:,col] = residue( seq )
    else
      throw("There is and aligned sequence with different number of columns [ $(length(seq.data)) != $(nres) ]:\n$(ascii(seq))")
    end
  end
  aln'
end

function __to_msa_mapping(sequences::Array{ASCIIString,1})
  nseq = size(sequences,1)
  nres = length(sequences[1])
  aln = Array(Residue,nres,nseq)
  mapp = Array(Int,nseq,nres)
  gaps = Uint8['.', '-']
  for i in 1:nseq
    init = 1
    seq = sequences[i].data
    if length(seq) == nres
      @inbounds for j in 1:nres
        res = seq[j]
        aln[j,i] = residue( res )
        if !( res in gaps )
          mapp[i, j] = init
          init += 1
        end
      end
    else
      throw("There is and aligned sequence with different number of columns [ $(length(seq.data)) != $(nres) ]:\n$(ascii(seq))")
    end
  end
  (aln', mapp)
end

function __to_msa_mapping(sequences::Array{ASCIIString,1}, ids::Array{ASCIIString,1})
  nseq = size(sequences,1)
  nres = length(sequences[1])
  aln = Array(Residue,nres,nseq)
  mapp = Array(Int,nseq,nres)
  gaps = Uint8['.', '-']
  sep = r"/|-"
  for i in 1:nseq
    fields = split(ids[i],sep)
    init = length(fields) == 3 ? int(fields[2]) : 1
    seq = sequences[i].data
    if length(seq) == nres
      @inbounds for j in 1:nres
        res = seq[j]
        aln[j,i] = residue( res )
        if !( res in gaps )
          mapp[i, j] = init
          init += 1
        end
      end
    else
      throw("There is and aligned sequence with different number of columns [ $(length(seq.data)) != $(nres) ]:\n$(ascii(seq))")
    end
    if (init - 1) != int(fields[3])
      throw("Different lengths: $(fields[3]) != $(init) for sequence $(ids[i])")
    end
  end
  (aln', mapp)
end

function readpfam(filename::ASCIIString; useidcoordinates::Bool=true)
  IDS, SEQS, GF, GS, GC, GR = __pre_readstockholm(filename)
  MSA, MAP = useidcoordinates ? __to_msa_mapping(SEQS, IDS) : __to_msa_mapping(SEQS)
  COLS = vcat(1:size(MSA,2))
  msa = MultipleSequenceAlignment(indexedvector(IDS), MSA, MAP, indexedvector(COLS), Annotations(GF, GS, GC, GR))
  filtercolumns!(msa, columngappercentage(msa) .< 1.0)
end

# Print Pfam
# ==========

function __printfileannotations(io::IO, msa::MultipleSequenceAlignment)
	if !isempty(msa.annotations.file)
		for (key, value) in msa.annotations.file
			println(io, string("#=GF ", key, " ", value))
		end
	end
end

function __printcolumnsannotations(io::IO, msa::MultipleSequenceAlignment)
	if !isempty(msa.annotations.columns)
		for (key, value) in msa.annotations.columns
			println(io, string("#=GC ", key, "\t\t", value))
		end
	end
end

function __printsequencesannotations(io::IO, msa::MultipleSequenceAlignment)
	if !isempty(msa.annotations.sequences)
		for (key, value) in msa.annotations.sequences
			println(io, string("#=GC ", key[1], " ", key[2], " ", value))
		end
	end
end


function __to_sequence_dict(annotation::Dict{(ASCIIString,ASCIIString),ASCIIString})
	seq_dict = Dict{ASCIIString,Vector{ASCIIString}}()
	for (key, value) in annotation
		seq_id = key[1]
		if seq_id in keys(seq_dict)
			push!(seq_dict[seq_id], string(key[2], "\t", value))
		else
			seq_dict[seq_id] = [ string(key[2], "\t", value) ]
		end
	end
	sizehint(seq_dict, length(seq_dict))
end

function printpfam(io::IO, msa::MultipleSequenceAlignment)
	__printfileannotations(io, msa)
	__printsequencesannotations(io, msa)
	res_annotations = __to_sequence_dict(msa.annotations.residues)
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
	__printcolumnsannotations(io, msa)
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
    throw(string(pfamcode, " is not a correct Pfam code"))
  end
end
