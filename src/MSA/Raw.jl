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
