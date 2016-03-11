# sub(aln, :, i) is ~ 0.015 seconds faster than aln[:,i] and better with the memmory (PF00085)
# aln.msa' takes ~ 0.015 seconds in PF00085
# Using aln.msa instead of aln in estimate_on_column_pairs is ~ 0.53 seconds faster for PF00085

"""
`estimateincolumns(aln, [count,] use, [α, β,] measure, [pseudocount,] [weight,] [usediagonal, diagonalvalue])`

This function `estimate` a `AbstractMeasure` in columns or pair of columns of a MSA.

- `aln` : This argument is mandatory and it should be a `Matrix{Residue}`. Use the function `getresidues` (from the MSA module) over a MSA object to get the needed matrix.
- `count` : This argument is optional. It should be defined when `use` is a `ResidueProbability` object. It indicates the element type of the counting table.
- `use` : This argument is mandatory and indicates the sub-type of `ResidueContingencyTables` used by `estimate` inside the function.
If the table has one dimension (`N`=`1`), the occurrences/probabilities are counted for each sequence/column.
If the table has two dimension (`N`=`2`), pairs of sequences/columns are used.
The dimension `N` and the `UseGap` parameter of `Residueount{T, N, UseGap}` or `ResidueProbability{T, N, UseGap}` determines the output and behaviour of this functions.
If `UseGap` is true, gaps are used in the estimations.
- `α` : This argument is optional, and indicates the weight of real frequencies to apply BLOSUM62 based pseudo frequencies.
- `β` : This argument is optional, and indicates the weight of BLOSUM62 based pseudo frequencies.
- `measure` : This argument is mandatory and indicates the measure to be used by `estimate` inside the function.
- `pseudocount` : This argument is optional. It should be an `AdditiveSmoothing` instance (default to zero).
- `weight` : This argument is optional. It should be an instance of `ClusteringResult` or `AbstractVector` (vector of weights).
Each sequence has weight 1 (`NoClustering()`) by default.
- `usediagonal` : This functions return a `Vector` in the one dimensional case, or a `PairwiseListMatrix` in the bidimensional case.
This argument only have sense in the bidimensional case and indicates if the list on the `PairwiseListMatrix` should include the diagonal (default to `true`).
- `diagonalvalue` : This argument is optional (default to zero). Indicates the value of output diagonal elements.
"""
function estimateincolumns{T, TP, UseGap}(aln::Matrix{Residue}, use::Type{ResidueCount{T, 1, UseGap}}, measure::AbstractMeasure{TP},
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

function estimateincolumns{T <: Real, TP, UseGap}(aln::Matrix{Residue}, count::Type{T}, use::Type{ResidueProbability{TP, 1, UseGap}}, measure::AbstractMeasure{TP},
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


"""
This function `estimate` a `measure` over sequences or sequence pairs.
It has the same arguments than `estimateincolumns`, look the documentation of the last.
"""
estimateinsequences(aln::Matrix{Residue}, args...) = estimateincolumns(transpose(aln), args...)
estimateinsequences(aln::AbstractMultipleSequenceAlignment, args...) = estimateincolumns(transpose(aln.msa), args...)
