# Calculates percent identity of two aligned sequences
# No account of the positions with gaps in both sequences in the length
function percentidentity(seq1, seq2)
  len = length(seq1)
  if len != length(seq2)
     throw("Sequences of different length, they aren't aligned or don't come from the same alignment")
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
