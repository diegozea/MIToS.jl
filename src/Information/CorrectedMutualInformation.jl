# function calculatezscore{T,N}(value::AbstractArray{T,N}, average::AbstractArray{T,N}, sd::AbstractArray{T,N}, unit::T=one(T)) # tolerance::T=eps(T))
#   zscore = similar(value)
#   if size(value) == size(average) == size(sd)
#     for i in eachindex(zscore)
#       val = value[i]
#       ave = average[i]
#       sta = sd[i]
#       if isapprox(val, ave) # abs(val - ave)  <= tolerance
#         zscore[i] = 0.0
#       elseif !isapprox(sta + unit, unit) # abs(sta) > tolerance
#         # Test for 0.0 using 1.0: 0.0 + 1.0 == 1.0
#         zscore[i] = (val - ave)/sta
#       else
#         zscore[i] = NaN
#       end
#     end
#   else
#     throw(ErrorException("The elements should have the same size"))
#   end
#   zscore
# end

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
  - usegap      Bool      false   Exclude gaps from statistics
  - samples     Int       100     Number of samples for Z-score
  - fixedgaps   Bool      true    Fix gaps positions for the random samples


This function returns 5 arrays:

  - Z score
  - MI or MIp
  - Mean of MI/MIp on the random samples
  - Standard deviation of MI/MIp on the random samples
  - List of used columns
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
