# Busjle et. al. 2009
# ===================

function _buslje09(aln, usegap, clusters, lambda, apc)
  mi = estimateincolumns(aln, ResidueCount{Float64, 2, usegap}, MutualInformation{Float64}(),
                         AdditiveSmoothing{Float64}(lambda), clusters, false)
  if apc
    APC!(mi)
  end
  mi
end

function _buslje09(aln::Matrix{Residue}; lambda::Float64=0.05,
                               clustering::Bool=true, threshold::Float64=0.62,
                               maxgap::Float64=0.5, apc::Bool=true, samples::Int=100,
                               usegap::Bool=false, fixedgaps::Bool=true)
  used = gappercentage(aln,1) .<= maxgap
  ncol = ncolumns(aln)
  aln = filtercolumns(aln, used)
  clusters = clustering ? hobohmI(aln, threshold) : NoClustering()
  mi = _buslje09(aln, usegap, clusters, lambda, apc)
  #rand_mi = Array(Float64, size(mi, 1), size(mi, 2), samples)
  rand_mi = Array(typeof(mi), samples)
  for ns in 1:samples
    fixedgaps ? shuffle_residues_sequencewise!(aln) : shuffle_sequencewise!(aln)
    rand_mi[ns] = _buslje09(aln, usegap, clusters, lambda, apc)
  end
  usedcol = collect(1:ncol)[used]
  (zscore(rand_mi, mi), mi, usedcol)
end

"""
This function takes a MSA or a file and a `Format` as first arguments.
Calculates a Z score and a corrected MI/MIp as described on **Busjle et. al. 2009**

Argument, type, default value and descriptions:

  - lambda      Float64   0.05    Low count value
  - clustering  Bool      true    Sequence clustering (Hobohm I)
  - threshold   Float64   0.62    Percent identity threshold for clustering
  - maxgap      Float64   0.5     Maximum fraction of gaps in positions included in calculation
  - apc         Bool      true    Use APC correction (MIp)
  - usegap      Bool      false   Use gaps on statistics
  - samples     Int       100     Number of samples for Z-score
  - fixedgaps   Bool      true    Fix gaps positions for the random samples


This function returns:

  - Z score
  - MI or MIp
"""
function buslje09(aln::Matrix{Residue}; kargs...)
  zscore, mi, used = _buslje09(aln; kargs...)
  labels!(mi, used)
  labels!(zscore, used)
  (zscore, mi)
end

buslje09(aln::MultipleSequenceAlignment; kargs...) = buslje09(aln.msa; kargs...)

function buslje09(aln::AnnotatedMultipleSequenceAlignment; kargs...)
  zscore, mi, used = _buslje09(aln.msa; kargs...)
  if haskey(getannotfile(aln), "ColMap")
    labels!(mi, getcolumnmapping(aln)[used])
    labels!(zscore, getcolumnmapping(aln)[used])
  else
    labels!(mi, used)
    labels!(zscore, used)
  end
  (zscore, mi)
end

function buslje09{T <: Format}(filename::AbstractString, format::Type{T}; kargs...)
  aln = read(filename, T, AnnotatedMultipleSequenceAlignment, generatemapping=true)
  buslje09(aln; kargs...)
end

# MIToS BLMI: Blosum MI
# =====================

function _BLMI(aln, clusters, alpha, beta, apc)
  mi = estimateincolumns(aln, Float64, ResidueProbability{Float64, 2, false}, alpha, beta, MutualInformation{Float64}(), zero(AdditiveSmoothing{Float64}) , clusters, false, 0.0)
  if apc
    APC!(mi)
  end
  mi
end

function _BLMI(aln::Matrix{Residue}; beta::Float64=4.6, threshold::Float64=0.62,
                                     maxgap::Float64=0.5, apc::Bool=true, samples::Int=100,
                                     fixedgaps::Bool=true)
  used = gappercentage(aln,1) .<= maxgap
  ncol = ncolumns(aln)
  aln = filtercolumns(aln, used)
  clusters = hobohmI(aln, threshold)
  numbercl = getnclusters(clusters)
  mi = _BLMI(aln, clusters, numbercl, beta, apc)
  rand_mi = Array(typeof(mi), samples)
  for ns in 1:samples
    fixedgaps ? shuffle_residues_sequencewise!(aln) : shuffle_sequencewise!(aln)
    rand_mi[ns] = _BLMI(aln, clusters, numbercl, beta, apc)
  end
  usedcol = collect(1:ncol)[used]
  (zscore(rand_mi, mi), mi, usedcol)
end

"""
This function takes a MSA or a file and a `Format` as first arguments.
Calculates a Z score (BLMI) and a corrected MI/MIp as described on **Busjle et. al. 2009** but using using BLOSUM62 pseudo frequencies instead of a fixed pseudocount.

Argument, type, default value and descriptions:

  - beta        Float64   4.6     β for BLOSUM62 pseudo frequencies
  - threshold   Float64   0.62    Percent identity threshold for sequence clustering (Hobohm I)
  - maxgap      Float64   0.5     Maximum fraction of gaps in positions included in calculation
  - apc         Bool      true    Use APC correction (MIp)
  - samples     Int       100     Number of samples for Z-score
  - fixedgaps   Bool      true    Fix gaps positions for the random samples


This function returns:

  - Z score (BLMI)
  - MI or MIp using BLOSUM62 pseudo frequencies
"""
function BLMI(aln::Matrix{Residue}; kargs...)
  zscore, mi, used = _BLMI(aln; kargs...)
  labels!(mi, used)
  labels!(zscore, used)
  (zscore, mi)
end

BLMI(aln::MultipleSequenceAlignment; kargs...) = BLMI(aln.msa; kargs...)

function BLMI(aln::AnnotatedMultipleSequenceAlignment; kargs...)
  zscore, mi, used = _BLMI(aln.msa; kargs...)
  if haskey(getannotfile(aln), "ColMap")
    labels!(mi, getcolumnmapping(aln)[used])
    labels!(zscore, getcolumnmapping(aln)[used])
  else
    labels!(mi, used)
    labels!(zscore, used)
  end
  (zscore, mi)
end

function BLMI{T <: Format}(filename::AbstractString, format::Type{T}; kargs...)
  aln = read(filename, T, AnnotatedMultipleSequenceAlignment, generatemapping=true)
  BLMI(aln; kargs...)
end

# MIToS Pairwise Gap Percentage
# =============================

function _pairwisegappercentage(aln::Matrix{Residue}; clustering::Bool=true, threshold::Float64=0.62)
  ncol = ncolumns(aln)
  clusters = clustering ? hobohmI(aln, threshold) : NoClustering()
  gu = estimateincolumns(aln, ResidueCount{Float64, 2, true}, GapUnionPercentage{Float64}(), zero(AdditiveSmoothing{Float64}), clusters, true)
  gi = estimateincolumns(aln, ResidueCount{Float64, 2, true}, GapIntersectionPercentage{Float64}(), zero(AdditiveSmoothing{Float64}), clusters, true)
  usedcol = collect(1:ncol)
  (gu, gi, usedcol)
end

"""
This function takes a MSA or a file and a `Format` as first arguments.
Calculates the percentage of gaps on columns pairs (union and intersection) using sequence clustering (Hobohm I).

Argument, type, default value and descriptions:

  - clustering  Bool      true    Sequence clustering (Hobohm I)
  - threshold   Float64   0.62    Percent identity threshold for sequence clustering (Hobohm I)

This function returns:

  - pairwise gap percentage (union)
  - pairwise gap percentage (intersection)
"""
function pairwisegappercentage(aln::Matrix{Residue}; kargs...)
  gu, gi, used = _pairwisegappercentage(aln; kargs...)
  labels!(gu, used)
  labels!(gi, used)
  (gu, gi)
end

pairwisegappercentage(aln::MultipleSequenceAlignment; kargs...) = pairwisegappercentage(aln.msa; kargs...)

function pairwisegappercentage(aln::AnnotatedMultipleSequenceAlignment; kargs...)
  gu, gi, used = _pairwisegappercentage(aln.msa; kargs...)
  if haskey(getannotfile(aln), "ColMap")
    map = getcolumnmapping(aln)
    labels!(gu, map)
    labels!(gi, map)
  else
    labels!(gu, used)
    labels!(gi, used)
  end
  (gu, gi)
end

function pairwisegappercentage{T <: Format}(filename::AbstractString, format::Type{T}; kargs...)
  aln = read(filename, T, AnnotatedMultipleSequenceAlignment, generatemapping=true)
  pairwisegappercentage(aln; kargs...)
end





