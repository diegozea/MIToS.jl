function _buslje09(aln, usegap, clusters, lambda, apc)
  mi = estimateincolumns(aln, ResidueCount{Float64, 2, usegap}, MutualInformation{Float64}(),
                         AdditiveSmoothing{Float64}(lambda), clusters)
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
  aln = filtercolumns(aln, gappercentage(aln,1) .<= maxgap)
  clusters = clustering ? hobohmI(aln, threshold) : 1
  mi = _buslje09(aln, usegap, clusters, lambda, apc)
  rand_mi = Array(Float64, size(mi, 1), size(mi, 2), samples)
  for ns in 1:samples
    fixedgaps ? shuffle_residues_sequencewise!(aln) : shuffle_sequencewise!(aln)
    rand_mi[:,:,ns] = _buslje09(aln, usegap, clusters, lambda, apc)
  end
  rand_mean = squeeze(mean(rand_mi,3),3)
  rand_sd = squeeze(std(rand_mi,3),3)
  (mi, (mi .- rand_mean) ./ rand_sd)
end

buslje09(aln::AnnotatedMultipleSequenceAlignment; kargs...) = buslje09(aln.msa; kargs...)

function buslje09{T <: Format}(filename::AbstractString, format::Type{T}; kargs...)
  aln = read(filename, T, AnnotatedMultipleSequenceAlignment, generatemapping=true)
  mi, zscore = buslje09(aln; kargs...)
  (getcolumnmapping(aln), mi, zscore)
end
