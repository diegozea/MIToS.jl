# Calculates percent identity of two aligned sequences
# No account of the positions with gaps in both sequences in the length
"""
seq1 and seq2 should have the same len
"""
function _percentidentity(seq1, seq2, len)
  count = zero(Int)
  colgap = zero(Int)
  @inbounds for i in 1:len
    if seq1[i] == seq2[i]
      count += one(Int)
      colgap += Int(seq1[i] == GAP)
    end
  end
  (count-colgap)/(len-colgap)
end

"""
Calculates the identity value between two aligned sequences.

The identity value is calculated as the number of identical characters in the i-th position of both
sequences divided by the length of both sequences.
Positions with gaps in both sequences are not counted in the length of the sequence.
Returns a value in [0, 1] range.
"""
function percentidentity(seq1, seq2)
  len = length(seq1)
  if len != length(seq2)
     throw(ErrorException("Sequences of different length, they aren't aligned or don't come from the same alignment"))
  end
  _percentidentity(seq1, seq2, len)
end

"""
Computes quickly if two aligned sequences have a identity value greater than a given `threshold` value.
`threshold` should be a number in [0, 1] range.
Returns a boolean value.
"""
function percentidentity(seq1, seq2, threshold::Float64)
  len = length(seq1)
  if len != length(seq2)
     throw("Sequences of different length, they aren't aligned or don't come from the same alignment")
  end
  n = len
  limit_count = n * threshold
  diff = 0
  count = 0
  for i in 1:len
    if seq1[i] == seq2[i]
      if seq1[i] != GAP
        count += 1
        if count >= limit_count
          return(true)
        end
      else
        n -= 1
        limit_count = n * threshold
      end
    else
      diff += 1
      if diff > n - limit_count
        return(false)
      end
    end
  end
  (count/n) >= threshold
end

# percentidentity for a MSA
# =========================

"""
aln should be transpose(msa)
"""
function _percentidentity_kernel!(scores, aln, nseq, len)
  k = 0
  list = scores.list
  @inbounds for i in 1:(nseq-1)
    a = aln[i]
    for j in (i+1):nseq
      list[k += 1] = _percentidentity(a, aln[j], len) * 100.0
    end
  end
  scores
end

"""
Calculates the identity between all the sequences on a MSA.
You can indicate the output element type with the last optional parameter (`Float64` by default).
For a MSA with a lot of sequences, you can use `Float32` or `Flot16` in order to avoid the `OutOfMemoryError()`.
"""
function percentidentity{T}(msa::AbstractMatrix{Residue}, out::Type{T}=Float64)
  aln = getresiduesequences(msa)
  nseq = length(aln)
  len = length(aln[1])
  scores = PairwiseListMatrix(T, nseq, false, T(100.0))
  _percentidentity_kernel!(scores, aln, nseq, len)
  scores
end
