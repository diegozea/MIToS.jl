function __deleteitems!(vector::Vector, items)
  i = 1
  j = 0
  while j < length(vector)
    j += 1
    if !(vector[j] in items)
      if i != j
        vector[i] = vector[j]
      end
      i += 1
    end
  end
  resize!(vector, i-1)
end

function __pre_readfasta(filename::ASCIIString)
  IDS  = ASCIIString[]
  SEQS = ASCIIString[]

  delim = UInt8('>')
  newline = UInt8('\n')
  deletechars = IntSet(Int[delim, newline, '\t', ' '])

  fh = open(filename,"r")
  seqinfo = readuntil(fh, delim)
  if seqinfo[1] == delim
    seqinfo = readuntil(fh, delim)
    while length(seqinfo) != 0
      firstnewline = findfirst(seqinfo, newline)
      push!(IDS, ascii( seqinfo[1:(firstnewline-1)] ) )
      push!(SEQS, ascii( __deleteitems!(seqinfo[(firstnewline+1):end], deletechars) ) )
      seqinfo = readuntil(fh, delim)
    end
  else
    throw("The file doesn't start with '>'")
  end
  close(fh)

  (IDS, SEQS)
end

function readfasta(filename::ASCIIString; useidcoordinates::Bool=true)
  IDS, SEQS = __pre_readfasta(filename)
  MSA, MAP = useidcoordinates ? __to_msa_mapping(SEQS, IDS) : __to_msa_mapping(SEQS)
  COLS = vcat(1:size(MSA,2))
  msa = MultipleSequenceAlignment(indexedvector(IDS), MSA, MAP, indexedvector(COLS), __empty(Annotations))
  filtercolumns!(msa, columngappercentage(msa) .< 1.0)
end

# Print Pfam
# ==========

function printfasta(io::IO, msa::MultipleSequenceAlignment)
	for i in 1:nsequences(msa)
		id = selectvalue(msa.id, i)
		seq = asciisequence(msa, i)
		println(io, string(">", id, "\n", seq))
	end
end

printfasta(msa::MultipleSequenceAlignment) = printfasta(STDOUT, msa)

# Write Pfam
# ==========

function writefasta(filename::ASCIIString, msa::MultipleSequenceAlignment)
	open(filename, "w") do fh
		printfasta(fh, msa)
	end
end

