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
      colgap += int(seq1[i] == GAP)
    end
  end
  (count-colgap)/(len-colgap)
end
