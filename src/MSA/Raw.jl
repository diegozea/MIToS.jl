immutable Raw <: Format end

# Raw Parser
# ==========

function parse(io::Union{IO, AbstractString}, format::Type{Raw}, output::Type{Matrix{Residue}}; deletefullgaps::Bool=true)
  SEQS = ASCIIString[]

  for line in eachline(io)
    push!(SEQS, chomp(line))
  end

  msa = convert(Matrix{Residue}, SEQS)

  if deletefullgaps
    return(deletefullgapcolumns(msa))
  end
  msa
end

function parse(io::Union{IO,AbstractString}, format::Type{Raw}, output::Type{AnnotatedMultipleSequenceAlignment}; generatemapping::Bool=false, deletefullgaps::Bool=true)
  SEQS = ASCIIString[]
  for line in eachline(io)
    push!(SEQS, chomp(line))
  end
  annot = Annotations()
  if generatemapping
    MSA, MAP = _to_msa_mapping(SEQS)
    setannotfile!(annot, "ColMap", join(vcat(1:size(MSA,2)), ','))
    for i in 1:length(SEQS)
      setannotsequence!(annot, string(i), "SeqMap", MAP[i])
    end
  else
    MSA = convert(Matrix{Residue}, SEQS)
  end
  msa = AnnotatedMultipleSequenceAlignment(IndexedArray(ASCIIString[ string(i) for i in 1:length(SEQS)]), MSA, annot)
  if deletefullgaps
    deletefullgapcolumns!(msa)
  end
  msa
end

parse(io::Union{IO, AbstractString}, format::Type{Raw}; deletefullgaps::Bool=true) = parse(io, Raw, Matrix{Residue}; deletefullgaps=deletefullgaps)

# Print Raw
# =========

function print(io::IO, msa::Matrix{Residue}, format::Type{Raw})
	for i in 1:nsequences(msa)
		seq = asciisequence(msa, i)
		println(io, seq)
	end
end

print(io::IO, msa::AbstractMultipleSequenceAlignment, format::Type{Raw}) = print(io, msa.msa, Raw)

print(msa::Matrix{Residue}, format::Type{Raw}) = print(STDOUT, msa, Raw)
print(msa::AbstractMultipleSequenceAlignment, format::Type{Raw}) = print(STDOUT, msa.msa, Raw)
