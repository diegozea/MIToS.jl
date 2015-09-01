immutable Count{T <: Real} end
immutable Probability{T <: Real} end

# sub(aln, :, i) is ~ 0.015 seconds faster than aln[:,i] and better with the memmory (PF00085)
# aln.msa' takes ~ 0.015 seconds in PF00085
# Using aln.msa instead of aln in estimate_on_column_pairs is ~ 0.53 seconds faster for PF00085

function estimate_on_columns{T}(aln::Matrix{Residue}, use::Type{Count{T}}, measure::InformationMeasure;
                                usegap::Bool=false, weight::SequenceWeights=1, pseudocount::Pseudocount=AdditiveSmoothing(zero(T)))
  N = ResidueCount{T, 1, usegap}()
  ncol = ncolumns(aln)
  scores = Array(Float64, ncol)
  for i in 1:ncol
    fill!(N, pseudocount) # instead of fill! with 0
    count!(N, weight, sub(aln,:,i))
    scores[i] = estimate(measure, N)
  end
  scores
end

function estimate_on_columns{T}(aln::Matrix{Residue}, use::Type{Probability{T}}, measure::InformationMeasure;
                                usegap::Bool=false, weight::SequenceWeights=1, pseudocount::Pseudocount=AdditiveSmoothing(zero(T)))
  N = ResidueCount{T, 1, usegap}()
  P = ResidueProbability{1, usegap}()
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

estimate_on_columns{T}(aln::AbstractMultipleSequenceAlignment, use::Type{Count{T}}, measure::InformationMeasure; usegap::Bool=false, weight::SequenceWeights=1,
                    pseudocount::Pseudocount=AdditiveSmoothing(zero(T))) = estimate_on_columns(aln.msa, use, measure; usegap=usegap, weight=weight, pseudocount=pseudocount)

estimate_on_columns{T}(aln::AbstractMultipleSequenceAlignment, use::Type{Probability{T}}, measure::InformationMeasure; usegap::Bool=false, weight::SequenceWeights=1,
                    pseudocount::Pseudocount=AdditiveSmoothing(zero(T))) = estimate_on_columns(aln.msa, use, measure; usegap=usegap, weight=weight, pseudocount=pseudocount)

function estimate_on_column_pairs{T}(aln::Matrix{Residue}, use::Type{Count{T}}, measure::SymmetricMeasure;
                                usegap::Bool=false, weight::SequenceWeights=1, pseudocount::Pseudocount=AdditiveSmoothing(zero(T)), diagonal::Float64=0.0)
  Nab = ResidueCount{T, 2, usegap}()
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

function estimate_on_column_pairs{T}(aln::Matrix{Residue}, use::Type{Probability{T}}, measure::SymmetricMeasure;
                                usegap::Bool=false, weight::SequenceWeights=1, pseudocount::Pseudocount=AdditiveSmoothing(zero(T)), diagonal::Float64=0.0)
  Nab = ResidueCount{T, 2, usegap}()
  Pab = ResidueProbability{2, usegap}()
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

estimate_on_column_pairs{T}(aln::AbstractMultipleSequenceAlignment, use::Type{Count{T}}, measure::InformationMeasure; usegap::Bool=false, weight::SequenceWeights=1,
                    pseudocount::Pseudocount=AdditiveSmoothing(zero(T)), diagonal::Float64=0.0) = estimate_on_columns(aln.msa, use, measure; usegap=usegap,
                                                                                                                      weight=weight, pseudocount=pseudocount, diagonal=diagonal)

estimate_on_column_pairs{T}(aln::AbstractMultipleSequenceAlignment, use::Type{Probability{T}}, measure::InformationMeasure; usegap::Bool=false, weight::SequenceWeights=1,
                    pseudocount::Pseudocount=AdditiveSmoothing(zero(T)), diagonal::Float64=0.0) = estimate_on_columns(aln.msa, use, measure; usegap=usegap,
                                                                                                                      weight=weight, pseudocount=pseudocount, diagonal=diagonal)

function estimate_on_column_pairs_blosum{T}(aln::Matrix{Residue}, use::Type{Probability{T}}, α, β, measure::SymmetricMeasure;
                                     weight::SequenceWeights=1, pseudocount::Pseudocount=AdditiveSmoothing(zero(T)), diagonal::Float64=0.0)
  Nab = ResidueCount{T, 2, false}()
  Pab = ResidueProbability{2, false}()
  Gab = ResidueProbability{2, false}()
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

function estimate_on_column_pairs_blosum{T}(aln::AbstractMultipleSequenceAlignment, use::Type{Probability{T}}, α, β, measure::SymmetricMeasure;
                                   weight::SequenceWeights=1, pseudocount::Pseudocount=AdditiveSmoothing(zero(T)), diagonal::Float64=0.0)
  estimate_on_column_pairs_blosum(aln.msa, use, α, β, measure; weight=weight, pseudocount=pseudocount, diagonal=diagonal)
end
