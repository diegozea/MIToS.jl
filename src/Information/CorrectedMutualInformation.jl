function calculatezscore{T,N}(value::AbstractArray{T,N}, average::AbstractArray{T,N}, sd::AbstractArray{T,N}, unit::T=one(T)) # tolerance::T=eps(T))
  zscore = similar(value)
  if size(value) == size(average) == size(sd)
    for i in eachindex(zscore)
      val = value[i]
      ave = average[i]
      sta = sd[i]
      if isapprox(val, ave) # abs(val - ave)  <= tolerance
        zscore[i] = 0.0
      elseif !isapprox(sta + unit, unit) # abs(sta) > tolerance
        # Test for 0.0 using 1.0: 0.0 + 1.0 == 1.0
        zscore[i] = (val - ave)/sta
      else
        zscore[i] = NaN
      end
    end
  else
    throw(ErrorException("The elements should have the same size"))
  end
  zscore
end

function _buslje09(aln, usegap, clusters, lambda, apc)
  mi = estimateincolumns(aln, ResidueCount{Float64, 2, usegap}, MutualInformation{Float64}(),
                         AdditiveSmoothing{Float64}(lambda), clusters, false)
  if apc
    APC!(mi)
  end
  mi
end

"""
#	[-c]                 1                    Sequence clustering [default on]
#	[-i float]           0.620000             Percent id for clustering
#	[-lc float]          0.050000             Low count value
#	[-ns int]            100                  Number of samples for Z-score

# [-xg]                1                    Exclude gaps from statistics [default on]
# [-fixg]              1                    Fixed gaps [default on]

# [-apc]               1                    Use APC correction [default on]
#	[-maxgap float]      0.500000             Max fraction of gaps in positions included in calculation
"""
function buslje09(aln::Matrix{Residue}; lambda::Float64=0.05,
                               clustering::Bool=true, threshold::Float64=0.62,
                               maxgap::Float64=0.5, apc::Bool=true, samples::Int=100,
                               usegap::Bool=false, fixedgaps::Bool=true)
  used = gappercentage(aln,1) .<= maxgap
  aln = filtercolumns(aln, used)
  clusters = clustering ? hobohmI(aln, threshold) : NoClustering()
  mi = _buslje09(aln, usegap, clusters, lambda, apc)
  rand_mi = Array(Float64, size(mi, 1), size(mi, 2), samples)
  for ns in 1:samples
    fixedgaps ? shuffle_residues_sequencewise!(aln) : shuffle_sequencewise!(aln)
    rand_mi[:,:,ns] = _buslje09(aln, usegap, clusters, lambda, apc)
  end
  rand_mean = squeeze(mean(rand_mi,3),3)
  rand_sd = squeeze(std(rand_mi,3),3)
  (collect(1:ncolumns(aln))[used],mi, rand_mean, rand_sd, calculatezscore(mi, rand_mean, rand_sd))
end

buslje09(aln::MultipleSequenceAlignment; kargs...) = buslje09(aln.msa; kargs...)

function buslje09(aln::AnnotatedMultipleSequenceAlignment; kargs...)
  used, mi, rand_mean, rand_sd, zscore = buslje09(aln.msa; kargs...)
  (getcolumnmapping(aln)[used], mi, rand_mean, rand_sd, zscore)
end

function buslje09{T <: Format}(filename::AbstractString, format::Type{T}; kargs...)
  aln = read(filename, T, AnnotatedMultipleSequenceAlignment, generatemapping=true)
  buslje09(aln; kargs...)
end
