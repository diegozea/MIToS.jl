# sub(aln, :, i) is ~ 0.015 seconds faster than aln[:,i] and better with the memmory (PF00085)
# aln.msa' takes ~ 0.015 seconds in PF00085
# Using aln.msa instead of aln in estimate_on_column_pairs is ~ 0.53 seconds faster for PF00085

function estimateincolumns{T, UseGap}(aln::Matrix{Residue}, use::Type{ResidueCount{T, 1, UseGap}}, measure::InformationMeasure;
                                      weight::SequenceWeights=1, pseudocount::Pseudocount=AdditiveSmoothing(zero(T)))
  N = ResidueCount{T, 1, UseGap}()
  ncol = ncolumns(aln)
  scores = Array(Float64, ncol)
  for i in 1:ncol
    fill!(N, pseudocount) # instead of fill! with 0
    count!(N, weight, sub(aln,:,i))
    scores[i] = estimate(measure, N)
  end
  scores
end

function estimateincolumns{T, TP, UseGap}(aln::Matrix{Residue}, count::Type{T}, use::Type{ResidueProbability{TP, 1, UseGap}}, measure::InformationMeasure;
                                usegap::Bool=false, weight::SequenceWeights=1, pseudocount::Pseudocount=AdditiveSmoothing(zero(T)))
  N = ResidueCount{T, 1, UseGap}()
  P = ResidueProbability{TP, 1, UseGap}()
  ncol = ncolumns(aln)
  scores = Array(Float64, ncol)
  for i in 1:ncol
    fill!(N, pseudocount) # instead of fill! with 0
    count!(N, weight, sub(aln,:,i))
    fill!(P, N) # count! calls update!
    scores[i] = estimate(measure, P)
  end
  scores
end

function estimateincolumns{T, UseGap}(aln::Matrix{Residue}, use::Type{ResidueCount{T, 2, UseGap}}, measure::SymmetricMeasure;
                                weight::SequenceWeights=1, pseudocount::Pseudocount=AdditiveSmoothing(zero(T)), diagonal::Float64=0.0)
  Nab = ResidueCount{T, 2, UseGap}()
  ncol = ncolumns(aln)
  scores = Array(Float64, ncol, ncol)
  @inbounds for i in 1:(ncol-1)
    a = sub(aln,:,i)
    for j in (i+1):ncol
      b = sub(aln,:,j)
      fill!(Nab, pseudocount) # instead of fill! with 0
      count!(Nab, weight, a, b)
      score = estimate(measure, Nab)
      scores[i,j] = score
      scores[j,i] = score
    end
  end
  @inbounds scores[diagind(scores)] = diagonal
  scores
end

function estimateincolumns{T, TP, UseGap}(aln::Matrix{Residue}, count::Type{T}, use::Type{ResidueProbability{TP, 2, UseGap}}, measure::SymmetricMeasure;
                                weight::SequenceWeights=1, pseudocount::Pseudocount=AdditiveSmoothing(zero(T)), diagonal::Float64=0.0)
  Nab = ResidueCount{T, 2, UseGap}()
  Pab = ResidueProbability{TP, 2, UseGap}()
  ncol = ncolumns(aln)
  scores = zeros(Float64, ncol, ncol)
  @inbounds for i in 1:(ncol-1)
    a = sub(aln,:,i)
    for j in (i+1):ncol
      b = sub(aln,:,j)
      fill!(Nab, pseudocount) # instead of fill! with 0
      count!(Nab, weight, a, b)
      fill!(Pab, Nab)
      score = estimate(measure, Pab)
      scores[i,j] = score
      scores[j,i] = score
    end
  end
  @inbounds scores[diagind(scores)] = diagonal
  scores
end

function estimateincolumns{T, TP}(aln::Matrix{Residue}, count::Type{T}, use::Type{ResidueProbability{TP, 2, false}}, α, β, measure::SymmetricMeasure;
                                     weight::SequenceWeights=1, pseudocount::Pseudocount=AdditiveSmoothing(zero(T)), diagonal::Float64=0.0)
  Nab = ResidueCount{T, 2, false}()
  Pab = ResidueProbability{TP, 2, false}()
  Gab = ResidueProbability{TP, 2, false}()
  ncol = ncolumns(aln)
  scores = zeros(Float64, ncol, ncol)
  @inbounds for i in 1:(ncol-1)
    a = sub(aln,:,i)
    for j in (i+1):ncol
      b = sub(aln,:,j)
      fill!(Nab, pseudocount) # instead of fill! with 0
      count!(Nab, weight, a, b)
      fill!(Pab, Nab)
      blosum_pseudofrequencies!(Gab, Pab)
      apply_pseudofrequencies!(Pab, Gab, α, β)
      score = estimate(measure, Pab)
      scores[i,j] = score
      scores[j,i] = score
    end
  end
  @inbounds scores[diagind(scores)] = diagonal
  scores
end

estimateincolumns(aln::AbstractMultipleSequenceAlignment, args...; kargs...) = estimateincolumns(aln.msa, args...; kargs...)

estimateinsequences(aln::Matrix{Residue}, args...; kargs...) = estimateincolumns(transpose(aln), args...; kargs...)
estimateinsequences(aln::AbstractMultipleSequenceAlignment, args...; kargs...) = estimateincolumns(transpose(aln.msa), args...; kargs...)
