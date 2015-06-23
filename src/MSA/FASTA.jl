using FastaIO

function readfasta(filename::ASCIIString; useidcoordinates::Bool=true)
  rawaln = FastaIO.readfasta("PF00085_full.fasta") # Make faster: push! sequences and/or reshape!
  nseq = length(rawaln)
  nres = length(rawaln[1][2])
  aln = Array(Residue,nres,nseq)
  mapp = Array(Int,nseq,nres)
  gaps = Uint8['.', '-']
  sep = r"/|-"
  ids = Array(ASCIIString,nseq)
  for i in 1:nseq
    id,seq = rawaln[i]
    ids[i] = id
    if useidcoordinates
      fields = split(id,sep)
      init = length(fields) == 3 ? int(fields[2]) : 1
    else
      init = 1
    end
    for j in 1:nres
      res = seq[j]
      aln[j,i] = residue(res)
      if !( res in gaps )
        mapp[i, j] = init
        init += 1
      end
    end
    if useidcoordinates && (init - 1) != int(fields[3])
      throw("Different lengths: $(fields[3]) != $(init) for sequence $(ids[i])")
    end
  end
  msa = MultipleSequenceAlignment(indexedvector(ids), aln', mapp, indexedvector([1:nres]), __empty(Annotations))
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

