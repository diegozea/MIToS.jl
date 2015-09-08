# Calculates percent identity of two aligned sequences
# No account of the positions with gaps in both sequences in the length
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
  count = zero(Int)
  colgap = zero(Int)
  for i in 1:len
    if seq1[i] == seq2[i]
      count += one(Int)
      colgap += Int(seq1[i] == GAP)
    end
  end
  (count-colgap)/(len-colgap)
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
