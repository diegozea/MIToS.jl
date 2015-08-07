using MIToS.Utils

function _pre_readfasta(filename::ASCIIString)
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
      push!(IDS, ASCIIString( seqinfo[1:(firstnewline-1)] ) )
      push!(SEQS, ASCIIString( deleteitems!(seqinfo[(firstnewline+1):end], deletechars) ) )
      seqinfo = readuntil(fh, delim)
    end
  else
    throw("The file doesn't start with '>'")
  end
  close(fh)

  (IDS, SEQS)
end

function readfasta(filename::ASCIIString; useidcoordinates::Bool=true)
  IDS, SEQS = _pre_readfasta(filename)
  MSA, MAP = useidcoordinates  && hascoordinates(IDS[1]) ? _to_msa_mapping(SEQS, IDS) : _to_msa_mapping(SEQS)
  COLS = vcat(1:size(MSA,2))
  msa = MultipleSequenceAlignment(IndexedVector(IDS), MSA, MAP, IndexedVector(COLS), Annotations())
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

