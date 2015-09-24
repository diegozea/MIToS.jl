# sub(aln, :, i) is ~ 0.015 seconds faster than aln[:,i] and better with the memmory (PF00085)
# aln.msa' takes ~ 0.015 seconds in PF00085
# Using aln.msa instead of aln in estimate_on_column_pairs is ~ 0.53 seconds faster for PF00085

function estimateincolumns{T, TP, UseGap}(aln::Matrix{Residue}, use::Type{ResidueCount{T, 1, UseGap}}, measure::InformationMeasure{TP},
                                          pseudocount::Pseudocount{T}=zero(AdditiveSmoothing{T}), weight::SequenceWeights=NoClustering())
  N = ResidueCount{T, 1, UseGap}()
  ncol = ncolumns(aln)
  scores = Array(TP, ncol)
  for i in 1:ncol
    fill!(N, pseudocount) # instead of fill! with 0
    count!(N, weight, sub(aln,:,i))
    scores[i] = estimate(measure, N)
  end
  scores
end

function estimateincolumns{T <: Real, TP, UseGap}(aln::Matrix{Residue}, count::Type{T}, use::Type{ResidueProbability{TP, 1, UseGap}}, measure::InformationMeasure{TP},
                                                  pseudocount::Pseudocount{T}=zero(AdditiveSmoothing{T}), weight::SequenceWeights=NoClustering())
  N = ResidueCount{T, 1, UseGap}()
  P = ResidueProbability{TP, 1, UseGap}()
  ncol = ncolumns(aln)
  scores = Array(TP, ncol)
  for i in 1:ncol
    fill!(N, pseudocount) # instead of fill! with 0
    count!(N, weight, sub(aln,:,i))
    fill!(P, N) # count! calls update!
    scores[i] = estimate(measure, P)
  end
  scores
end

# function estimateincolumns{T, TP, UseGap}(aln::Matrix{Residue}, use::Type{ResidueCount{T, 2, UseGap}}, measure::SymmetricMeasure{TP},
#                                           pseudocount::Pseudocount{T}=zero(AdditiveSmoothing{T}), weight::SequenceWeights=NoClustering(), usediagonal::Bool=true, diagonalvalue::TP=zero(TP))
#   Nab = ResidueCount{T, 2, UseGap}()
#   ncol = ncolumns(aln)
#   scores = Array(TP, ncol, ncol)
#   @inbounds for i in 1:ncol
#     a = sub(aln,:,i)
#     for j in i:ncol
#       if !usediagonal && i == j
#         scores[i,j] = diagonalvalue
#         continue
#       end
#       b = sub(aln,:,j)
#       fill!(Nab, pseudocount) # instead of fill! with 0
#       count!(Nab, weight, a, b)
#       score = estimate(measure, Nab)
#       scores[i,j] = score
#       scores[j,i] = score
#     end
#   end
#   scores
# end

# function estimateincolumns{T <: Real, TP, UseGap}(aln::Matrix{Residue}, count::Type{T}, use::Type{ResidueProbability{TP, 2, UseGap}},
#                                                   measure::SymmetricMeasure{TP}, pseudocount::Pseudocount{T}=zero(AdditiveSmoothing{T}),
#                                                   weight::SequenceWeights=NoClustering(), usediagonal::Bool=true, diagonalvalue::TP=zero(TP))
#   Nab = ResidueCount{T, 2, UseGap}()
#   Pab = ResidueProbability{TP, 2, UseGap}()
#   ncol = ncolumns(aln)
#   scores = zeros(TP, ncol, ncol)
#   @inbounds for i in 1:ncol
#     a = sub(aln,:,i)
#     for j in i:ncol
#       if !usediagonal && i == j
#         scores[i,j] = diagonalvalue
#         continue
#       end
#       b = sub(aln,:,j)
#       fill!(Nab, pseudocount) # instead of fill! with 0
#       count!(Nab, weight, a, b)
#       fill!(Pab, Nab)
#       score = estimate(measure, Pab)
#       scores[i,j] = score
#       scores[j,i] = score
#     end
#   end
#   scores
# end

# function estimateincolumns{T <: Real, TP}(aln::Matrix{Residue}, count::Type{T}, use::Type{ResidueProbability{TP, 2, false}}, α, β,
#                                           measure::SymmetricMeasure{TP}, pseudocount::Pseudocount{T}=zero(AdditiveSmoothing{T}),
#                                           weight::SequenceWeights=NoClustering(), usediagonal::Bool=true, diagonalvalue::TP=zero(TP))
#   Nab = ResidueCount{T, 2, false}()
#   Pab = ResidueProbability{TP, 2, false}()
#   Gab = ResidueProbability{TP, 2, false}()
#   ncol = ncolumns(aln)
#   scores = zeros(TP, ncol, ncol)
#   @inbounds for i in 1:ncol
#     a = sub(aln,:,i)
#     for j in i:ncol
#       if !usediagonal && i == j
#         scores[i,j] = diagonalvalue
#         continue
#       end
#       b = sub(aln,:,j)
#       fill!(Nab, pseudocount) # instead of fill! with 0
#       count!(Nab, weight, a, b)
#       fill!(Pab, Nab)
#       blosum_pseudofrequencies!(Gab, Pab)
#       apply_pseudofrequencies!(Pab, Gab, α, β)
#       score = estimate(measure, Pab)
#       scores[i,j] = score
#       scores[j,i] = score
#     end
#   end
#   scores
# end

using PairwiseListMatrices

function estimateincolumns{T, TP, UseGap}(aln::Matrix{Residue}, use::Type{ResidueCount{T, 2, UseGap}}, measure::SymmetricMeasure{TP},
                                          pseudocount::Pseudocount{T}=zero(AdditiveSmoothing{T}), weight::SequenceWeights=NoClustering(), usediagonal::Bool=true, diagonalvalue::TP=zero(TP))
  Nab = ResidueCount{T, 2, UseGap}()
  ncol = ncolumns(aln)
  scores = PairwiseListMatrix(TP, ncol, usediagonal, diagonalvalue) # Array(TP, ncol, ncol)
  @inbounds for i in 1:ncol
    a = sub(aln,:,i)
    for j in i:ncol
      if !usediagonal && i == j
        continue
      end
      b = sub(aln,:,j)
      fill!(Nab, pseudocount) # instead of fill! with 0
      count!(Nab, weight, a, b)
      scores[i,j] = estimate(measure, Nab)
    end
  end
  scores
end

function estimateincolumns{T <: Real, TP, UseGap}(aln::Matrix{Residue}, count::Type{T}, use::Type{ResidueProbability{TP, 2, UseGap}},
                                                  measure::SymmetricMeasure{TP}, pseudocount::Pseudocount{T}=zero(AdditiveSmoothing{T}),
                                                  weight::SequenceWeights=NoClustering(), usediagonal::Bool=true, diagonalvalue::TP=zero(TP))
  Nab = ResidueCount{T, 2, UseGap}()
  Pab = ResidueProbability{TP, 2, UseGap}()
  ncol = ncolumns(aln)
  scores = PairwiseListMatrix(TP, ncol, usediagonal, diagonalvalue) # zeros(TP, ncol, ncol)
  @inbounds for i in 1:ncol
    a = sub(aln,:,i)
    for j in i:ncol
      if !usediagonal && i == j
        continue
      end
      b = sub(aln,:,j)
      fill!(Nab, pseudocount) # instead of fill! with 0
      count!(Nab, weight, a, b)
      fill!(Pab, Nab)
      scores[i,j] = estimate(measure, Pab)
    end
  end
  scores
end

function estimateincolumns{T <: Real, TP}(aln::Matrix{Residue}, count::Type{T}, use::Type{ResidueProbability{TP, 2, false}}, α, β,
                                          measure::SymmetricMeasure{TP}, pseudocount::Pseudocount{T}=zero(AdditiveSmoothing{T}),
                                          weight::SequenceWeights=NoClustering(), usediagonal::Bool=true, diagonalvalue::TP=zero(TP))
  Nab = ResidueCount{T, 2, false}()
  Pab = ResidueProbability{TP, 2, false}()
  Gab = ResidueProbability{TP, 2, false}()
  ncol = ncolumns(aln)
  scores = PairwiseListMatrix(TP, ncol, usediagonal, diagonalvalue) # zeros(TP, ncol, ncol)
  @inbounds for i in 1:ncol
    a = sub(aln,:,i)
    for j in i:ncol
      if !usediagonal && i == j
        continue
      end
      b = sub(aln,:,j)
      fill!(Nab, pseudocount) # instead of fill! with 0
      count!(Nab, weight, a, b)
      fill!(Pab, Nab)
      blosum_pseudofrequencies!(Gab, Pab)
      apply_pseudofrequencies!(Pab, Gab, α, β)
      scores[i,j] = estimate(measure, Pab)
    end
  end
  scores
end

estimateincolumns(aln::AbstractMultipleSequenceAlignment, args...) = estimateincolumns(aln.msa, args...)

estimateinsequences(aln::Matrix{Residue}, args...) = estimateincolumns(transpose(aln), args...)
estimateinsequences(aln::AbstractMultipleSequenceAlignment, args...) = estimateincolumns(transpose(aln.msa), args...)
