using MIToS.Utils

import MIToS.Utils: reader

immutable FASTA <: Format end

function _pre_readfasta(io::IO)
  IDS  = ASCIIString[]
  SEQS = ASCIIString[]

  delim = UInt8('>')
  newline = UInt8('\n')
  deletechars = IntSet(Int[delim, newline, '\t', ' '])

  seqinfo = readuntil(io, delim)
  if seqinfo[1] == delim
    seqinfo = readuntil(io, delim)
    while length(seqinfo) != 0
      firstnewline = findfirst(seqinfo, newline)
      push!(IDS, ASCIIString( seqinfo[1:(firstnewline-1)] ) )
      push!(SEQS, ASCIIString( deleteitems!(seqinfo[(firstnewline+1):end], deletechars) ) )
      seqinfo = readuntil(io, delim)
    end
  else
    throw(ErrorException("The file doesn't start with '>'"))
  end

  (IDS, SEQS)
end

function reader(io::IO, format::Type{FASTA}, output::Type{AnnotatedMultipleSequenceAlignment}; useidcoordinates::Bool=true, deletefullgaps::Bool=true)
  IDS, SEQS = _pre_readfasta(io)
  MSA, MAP = useidcoordinates  && hascoordinates(IDS[1]) ? _to_msa_mapping(SEQS, IDS) : _to_msa_mapping(SEQS)
  COLS = vcat(1:size(MSA,2))
  msa = AnnotatedMultipleSequenceAlignment(IndexedVector(IDS), MSA, MAP, IndexedVector(COLS), empty(Annotations))
  if deletefullgaps
    deletefullgaps!(msa)
  end
  return(msa)
end

function reader(io::IO, format::Type{FASTA}, output::Type{MultipleSequenceAlignment}; deletefullgaps::Bool=true)
  IDS, SEQS = _pre_readfasta(io)
  msa = MultipleSequenceAlignment(IndexedVector(IDS), convert(Matrix{Residue}, SEQS))
  if deletefullgaps
    deletefullgaps!(msa)
  end
  return(msa)
end

reader(io::IO, format::Type{FASTA}; useidcoordinates::Bool=true, deletefullgaps::Bool=true) = reader(io, FASTA, AnnotatedMultipleSequenceAlignment; useidcoordinates=useidcoordinates, deletefullgaps=deletefullgaps)

# Print Pfam
# ==========

function printfasta(io::IO, msa::AbstractMultipleSequenceAlignment)
	for i in 1:nsequences(msa)
		id = selectvalue(msa.id, i)
		seq = asciisequence(msa, i)
		println(io, string(">", id, "\n", seq))
	end
end

printfasta(msa::AbstractMultipleSequenceAlignment) = printfasta(STDOUT, msa)

# Write Pfam
# ==========

function writefasta(filename::ASCIIString, msa::AbstractMultipleSequenceAlignment)
	open(filename, "w") do fh
		printfasta(fh, msa)
	end
end

