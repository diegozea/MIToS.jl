immutable Count{T <: Real} end
immutable Probability{T <: Real} end

function estimate_sequences{T}(aln, use::Type{Count{T}}, measure::InformationMeasure; usegap::Bool=false, weight::SequenceWeights=1, pseudocount::Pseudocount=AdditiveSmoothing(zero(T)))
  N = ResidueCount{T, 1, usegap}()
  len = nsequences(aln)
  scores = Array(Float64, len)
  for i in 1:len
    fill!(N, pseudocount) # instead of fill! with 0
    count!(N, weight, vec(aln[i,:])) # Test: subarray, slice or view
    scores[i] = estimate(measure, N)
  end
  scores
end

function estimate_sequences{T}(aln, use::Type{Probability{T}}, measure::InformationMeasure; usegap::Bool=false, weight::SequenceWeights=1, pseudocount::Pseudocount=AdditiveSmoothing(zero(T)))
  N = ResidueCount{T, 1, usegap}()
  P = ResidueProbability{1, usegap}()
  len = nsequences(aln)
  scores = Array(Float64, len)
  for i in 1:len
    fill!(N, pseudocount) # instead of fill! with 0
    count!(N, weight, vec(aln[i,:])) # Test: subarray, slice or view
    fill!(P, count) # count! calls update!
    scores[i] = estimate(measure, P)
  end
  scores
end


# """This method use BLOSUM62 based pseudofrequencies"""
# function probabilities(α, β, res1::AbstractVector{Residue}, res2::AbstractVector{Residue}; weight::SequenceWeights=1)
# 	Pab = fill!(ResidueProbability{2, false}(), count(res1, res2, usegap=false, weight=weight))
# 	Gab = blosum_pseudofrequencies!(ResidueProbability{2,false}(), Pab)
# 	apply_pseudofrequencies!(Pab, Gab, α, β)
# end
