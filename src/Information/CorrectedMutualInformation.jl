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

"""
This function takes a MSA or a file and a `Format` as first arguments.
Calculates a Z score and a corrected MI/MIp as described on **Busjle et. al. 2009**

Argument, type, default value and descriptions:

```
  - lambda      Float64   0.05    Low count value
  - clustering  Bool      true    Sequence clustering (Hobohm I)
  - threshold             62      Percent identity threshold for clustering
  - maxgap      Float64   0.5     Maximum fraction of gaps in positions included in calculation
  - apc         Bool      true    Use APC correction (MIp)
  - usegap      Bool      false   Use gaps on statistics
  - samples     Int       100     Number of samples for Z-score
  - fixedgaps   Bool      true    Fix gaps positions for the random samples
```

This function returns:

```
  - Z score
  - MI or MIp
```
"""
function buslje09(aln::AbstractMatrix{Residue}; lambda::Float64=0.05,
                  clustering::Bool=true, threshold=62,
                  maxgap::Float64=0.5, apc::Bool=true, samples::Int=100,
                  usegap::Bool=false, fixedgaps::Bool=true)
    used = gapfraction(aln,1) .<= maxgap
    col_labels = columnlabels(aln)[used]
    ncol = ncolumns(aln)
    aln = filtercolumns(aln, used)
    clusters = clustering ? hobohmI(aln, threshold) : NoClustering()
    mi = _buslje09(aln, usegap, clusters, lambda, apc)
    labels!(mi, col_labels)
    if samples > 0
        rand_mi = Array(typeof(mi), samples)
        for ns in 1:samples
            fixedgaps ? shuffle_residues_sequencewise!(aln) : shuffle_sequencewise!(aln)
            rand_mi[ns] = _buslje09(aln, usegap, clusters, lambda, apc)
        end
        zmi = zscore(rand_mi, mi)
        labels!(zmi, col_labels)
        return(zmi, mi)
    else
        return(zeros(mi), mi)
    end
end

function buslje09{T <: Format}(filename::AbstractString, format::Type{T}; kargs...)
    aln = read(filename, T, AnnotatedMultipleSequenceAlignment, generatemapping=true)
    buslje09(aln; kargs...)
end

# MIToS BLMI: Blosum MI
# =====================

function _BLMI(aln, clusters, alpha, beta, apc, lambda::Float64=zero(Float64))
    mi = estimateincolumns(aln, Float64, ResidueProbability{Float64, 2, false},
                           alpha, beta,
                           MutualInformation{Float64}(),
                           AdditiveSmoothing(lambda), clusters, false, 0.0)
    if apc
        APC!(mi)
    end
    mi
end

"""
This function takes a MSA or a file and a `Format` as first arguments.
Calculates a Z score (ZBLMI) and a corrected MI/MIp as described on **Busjle et. al. 2009** but using using BLOSUM62 pseudo frequencies instead of a fixed pseudocount.

Argument, type, default value and descriptions:

```
  - beta        Float64   8.512   β for BLOSUM62 pseudo frequencies
  - lambda      Float64   0.0     Low count value
  - threshold             62      Percent identity threshold for sequence clustering (Hobohm I)
  - maxgap      Float64   0.5     Maximum fraction of gaps in positions included in calculation
  - apc         Bool      true    Use APC correction (MIp)
  - samples     Int       50      Number of samples for Z-score
  - fixedgaps   Bool      true    Fix gaps positions for the random samples
```

This function returns:

```
  - Z score (ZBLMI)
  - MI or MIp using BLOSUM62 pseudo frequencies (BLMI/BLMIp)
```
"""
function BLMI(aln::AbstractMatrix{Residue}; beta::Float64=8.512, threshold=62,
              maxgap::Float64=0.5, apc::Bool=true, samples::Int=50,
              fixedgaps::Bool=true, lambda::Float64=zero(Float64))
    used = gapfraction(aln,1) .<= maxgap
    col_labels = columnlabels(aln)[used]
    ncol = ncolumns(aln)
    aln = filtercolumns(aln, used)
    clusters = hobohmI(aln, threshold)
    numbercl = nclusters(clusters)
    mi = _BLMI(aln, clusters, numbercl, beta, apc, lambda)
    labels!(mi, col_labels)
    if samples > 0
        rand_mi = Array(typeof(mi), samples)
        for ns in 1:samples
            fixedgaps ? shuffle_residues_sequencewise!(aln) : shuffle_sequencewise!(aln)
            rand_mi[ns] = _BLMI(aln, clusters, numbercl, beta, apc, lambda)
        end
        zmi = zscore(rand_mi, mi)
        labels!(zmi, col_labels)
        return(zmi, mi)
    else
        return(zeros(mi), mi)
    end
end

function BLMI{T <: Format}(filename::AbstractString, format::Type{T}; kargs...)
    aln = read(filename, T, AnnotatedMultipleSequenceAlignment, generatemapping=true)
    BLMI(aln; kargs...)
end

# MIToS Pairwise Gap Percentage
# =============================

"""
This function takes a MSA or a file and a `Format` as first arguments.
Calculates the percentage of gaps on columns pairs (union and intersection) using sequence clustering (Hobohm I).

Argument, type, default value and descriptions:

  - clustering  Bool      true    Sequence clustering (Hobohm I)
  - threshold             62      Percent identity threshold for sequence clustering (Hobohm I)

This function returns:

  - pairwise gap percentage (union)
  - pairwise gap percentage (intersection)
"""
function pairwisegapfraction(aln::AbstractMatrix{Residue}; clustering::Bool=true, threshold=62)
    col_labels = columnlabels(aln)
    clusters = clustering ? hobohmI(aln, threshold) : NoClustering()
    gu = estimateincolumns(aln, ResidueCount{Float64, 2, true}, GapUnionPercentage{Float64}(), zero(AdditiveSmoothing{Float64}), clusters, true)
    gi = estimateincolumns(aln, ResidueCount{Float64, 2, true}, GapIntersectionPercentage{Float64}(), zero(AdditiveSmoothing{Float64}), clusters, true)
    labels!(gu, col_labels)
    labels!(gi, col_labels)
    gu, gi
end

function pairwisegapfraction{T <: Format}(filename::AbstractString, format::Type{T}; kargs...)
    aln = read(filename, T, AnnotatedMultipleSequenceAlignment, generatemapping=true)
    pairwisegapfraction(aln; kargs...)
end
