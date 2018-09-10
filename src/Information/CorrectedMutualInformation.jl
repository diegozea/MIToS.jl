# Busjle et. al. 2009
# ===================

const _MI_MAT_TYPE = NamedArray{Float64,2,PairwiseListMatrix{Float64,false,Vector{Float64}},NTuple{2,OrderedDict{String,Int}}}

function _buslje09(aln, alphabet::A, clusters, lambda, apc) where A
    mi = mapcolpairfreq!(mutual_information, aln,
                         Counts{Float64,2,A}(ContingencyTable(Float64,Val{2},alphabet)),
                         Val{false}, # use diagonal
                         pseudocounts = AdditiveSmoothing{Float64}(lambda),
                         weights = clusters,
                         diagonalvalue = NaN)
    if apc
        APC!(mi)
    end
    mi
end

"""
`buslje09` takes a MSA or a file and a `FileFormat` as first arguments. It calculates a Z score
and a corrected MI/MIp as described on **Busjle et. al. 2009**.

keyword argument, type, default value and descriptions:

```
  - lambda      Float64   0.05    Low count value
  - clustering  Bool      true    Sequence clustering (Hobohm I)
  - threshold             62      Percent identity threshold for clustering
  - maxgap      Float64   0.5     Maximum fraction of gaps in positions included in calculation
  - apc         Bool      true    Use APC correction (MIp)
  - samples     Int       100     Number of samples for Z-score
  - fixedgaps   Bool      true    Fix gaps positions for the random samples
  - alphabet    ResidueAlphabet UngappedAlphabet()  Residue alphabet to be used
```

This function returns:

```
  - Z score
  - MI or MIp
```
"""
function buslje09(aln::AbstractMatrix{Residue};
                  lambda::Float64 = 0.05,
                  clustering::Bool = true,
                  threshold = 62,
                  maxgap::Float64 = 0.5,
                  apc::Bool = true,
                  samples::Int = 100,
                  alphabet::ResidueAlphabet = UngappedAlphabet(),
                  fixedgaps::Bool = true)::NTuple{2,_MI_MAT_TYPE}
    aln = filtercolumns(aln, gapfraction(aln,1) .<= maxgap)
    clusters = clustering ? hobohmI(aln, threshold) : NoClustering()
    mi = _buslje09(aln, alphabet, clusters, lambda, apc)
    if samples > 0
        rand_mi = Array{PairwiseListMatrix{Float64,false,Vector{Float64}}}(undef, samples)
        zmi = copy(mi)
        residuematrix = getresidues(aln)
        for ns in 1:samples
            shuffle!(residuematrix,1,fixedgaps)
            rand_mi[ns] = getarray(_buslje09(aln, alphabet, clusters, lambda, apc))
        end
        PairwiseListMatrices.zscore!(rand_mi, getarray(zmi))
        return(zmi, mi)
    else
        return(fill!(copy(mi),0.0), mi)
    end
end

function buslje09(filename::String, format::Type{T}; kargs...) where T <: FileFormat
    aln = read(filename, T, AnnotatedMultipleSequenceAlignment, generatemapping=true)
    buslje09(aln; kargs...)
end

# MIToS BLMI: Blosum MI
# =====================

function _BLMI(aln, clusters, alpha, beta, apc, lambda::Float64=0.0)
    mi = mapcolpairfreq!(mutual_information, aln,
                         Probabilities{Float64,2,UngappedAlphabet}(
                            ContingencyTable(Float64,Val{2},UngappedAlphabet())
                         ),
                         Val{false}, # use diagonal
                         pseudocounts = AdditiveSmoothing{Float64}(lambda),
                         weights = clusters,
                         diagonalvalue = NaN,
                         pseudofrequencies=BLOSUM_Pseudofrequencies(alpha,beta))
    if apc
        APC!(mi)
    end
    mi
end

"""
`BLMI` takes a MSA or a file and a `FileFormat` as first arguments. It calculates a Z score
(ZBLMI) and a corrected MI/MIp as described on **Busjle et. al. 2009** but using using
BLOSUM62 pseudo frequencies instead of a fixed pseudocount.

Keyword argument, type, default value and descriptions:

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
function BLMI(aln::AbstractMatrix{Residue};
              beta::Float64 = 8.512,
              threshold = 62,
              maxgap::Float64 = 0.5,
              apc::Bool = true,
              samples::Int = 50,
              fixedgaps::Bool = true,
              lambda::Float64 = 0.0)::NTuple{2,_MI_MAT_TYPE}
    aln = filtercolumns(aln, gapfraction(aln,1) .<= maxgap)
    clusters = hobohmI(aln, threshold)
    numbercl = Float64(nclusters(clusters))
    mi = _BLMI(aln, clusters, numbercl, beta, apc, lambda)
    if samples > 0
        rand_mi = Array{PairwiseListMatrix{Float64,false,Vector{Float64}}}(undef, samples)
        zmi = copy(mi)
        residuematrix = getresidues(aln)
        for ns in 1:samples
            shuffle!(residuematrix,1,fixedgaps)
            rand_mi[ns] = getarray(_BLMI(aln, clusters, numbercl, beta, apc, lambda))
        end
        PairwiseListMatrices.zscore(rand_mi, getarray(zmi))
        return(zmi, mi)
    else
        return(fill!(copy(mi),0.0), mi)
    end
end

function BLMI(filename::String, format::Type{T}; kargs...) where T <: FileFormat
    aln = read(filename, T, AnnotatedMultipleSequenceAlignment, generatemapping=true)
    BLMI(aln; kargs...)
end
