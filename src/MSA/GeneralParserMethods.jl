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
      throw( ErrorException("There is and aligned sequence with different number of columns [ $(length(seq)) != $(nres) ]:\n$(ascii(seq))") )
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
      throw( ErrorException("There is and aligned sequence with different number of columns [ $(length(seq)) != $(nres) ]:\n$(ascii(seq))") )
    end
    if (init - 1) != parse(Int, fields[3])
      throw( ErrorException("Different lengths: $(fields[3]) != $(init) for sequence $(ids[i])") )
    end
  end
  (aln', mapp)
end

function deletefullgaps!(msa::AbstractMultipleSequenceAlignment)
  mask = columngappercentage(msa) .!= one(Float64)
  if sum(mask) != zero(Float64)
    filtercolumns!(msa, mask)
  end
  msa
end