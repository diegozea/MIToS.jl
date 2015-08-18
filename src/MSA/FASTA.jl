using MIToS.Utils

import Base: parse

immutable FASTA <: Format end

# FASTA Parser
# ============

function _pre_readfasta(io::AbstractString)
  seqs = split(io, '>')
  N = length(seqs) - 1

  IDS  = Array(ASCIIString, N)
  SEQS = Array(ASCIIString, N)

  for i in 1:N
    fields = split(seqs[i+1], '\n')
    IDS[i] = replace(fields[1], r"\s+", "")
    SEQS[i] = replace(fields[2], r"\s+", "")
  end

  (IDS, SEQS)
end

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

function parse(io::Union(IO,AbstractString), format::Type{FASTA}, output::Type{AnnotatedMultipleSequenceAlignment}; useidcoordinates::Bool=true, deletefullgaps::Bool=true)
  IDS, SEQS = _pre_readfasta(io)
  MSA, MAP = useidcoordinates  && hascoordinates(IDS[1]) ? _to_msa_mapping(SEQS, IDS) : _to_msa_mapping(SEQS)
  COLS = vcat(1:size(MSA,2))
  msa = AnnotatedMultipleSequenceAlignment(IndexedVector(IDS), MSA, MAP, IndexedVector(COLS), empty(Annotations))
  if deletefullgaps
    deletefullgaps!(msa)
  end
  return(msa)
end

function parse(io::Union(IO,AbstractString), format::Type{FASTA}, output::Type{MultipleSequenceAlignment}; deletefullgaps::Bool=true)
  IDS, SEQS = _pre_readfasta(io)
  msa = MultipleSequenceAlignment(IndexedVector(IDS), convert(Matrix{Residue}, SEQS))
  if deletefullgaps
    deletefullgaps!(msa)
  end
  return(msa)
end

parse(io::Union(IO,AbstractString), format::Type{FASTA}; useidcoordinates::Bool=true, deletefullgaps::Bool=true) = parse(io, FASTA, AnnotatedMultipleSequenceAlignment; useidcoordinates=useidcoordinates, deletefullgaps=deletefullgaps)

# Print FASTA
# ===========

function print(io::IO, msa::AbstractMultipleSequenceAlignment, format::Type{FASTA})
	for i in 1:nsequences(msa)
		id = selectvalue(msa.id, i)
		seq = asciisequence(msa, i)
		println(io, string(">", id, "\n", seq))
	end
end

print(msa::MultipleSequenceAlignment, format::Type{FASTA}) = printfasta(STDOUT, msa, FASTA)
