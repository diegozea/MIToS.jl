## Pab
function fill_pab!(pab::Array{Float64,2},
                   seqi::AbstractVector{Residue}, seqj::AbstractVector{Residue},
                   weight::Vector{Float64}, pseudocount::Float64)
  fill!( pab, pseudocount )
  len = length(seqi)
#  if length(seqj) != len || length(weight) != len
#    throw("Different lengths")
#  end
  suma = sum(pab)
  for i in 1:len
    if seqi[i] != GAP  && seqj[i] != GAP
      @inbounds pab[ seqi[i] , seqj[i] ] += weight[i]
      @inbounds suma += weight[i]
    end
  end
  pab[:,:] /= suma
end

## Pseudo frequencies
function pab_pseudofrequencies!(pab::Array{Float64,2}, gab::Array{Float64,2},
                                nclusters::Int, beta::Float64)
  if beta == 0.0
	  return( pab )
  end
  suma = 0.0
  for a in 1:20, b in 1:20
    @inbounds gab[a,b] = 0.0
    for i in 1:20, j in 1:20
      if pab[i,j] != 0
	      # BLOSUM62_P_i_j[i,a] is p(a | i)
	      @inbounds val = pab[i,j] * BLOSUM62_P_i_j[i,a] * BLOSUM62_P_i_j[j,b]
	      @inbounds gab[a,b] += val
	      suma += val
      end
    end
  end
  gab[:,:] /= suma
  alpha = nclusters
  frac = 1 / ( alpha + beta )
  suma = 0.0
  for i in 1:20, j in 1:20
    @inbounds val = ( alpha * pab[i,j] + beta * gab[i,j] ) * frac
    @inbounds pab[i,j] = val
    suma += val
  end
  pab[:,:] /= suma
  pab
end
