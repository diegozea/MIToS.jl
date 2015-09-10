"""
Mean mutual information of column a (Dunn et. al. 2008).
Summation is over j=1 to N, j ≠ a. Total is N-1.
"""
_mean_column(mi::Matrix) = (squeeze(sum(mi, 1),1) .- diag(mi)) ./ (size(mi,1)-1)

"""
Mean mutual information of column a (Dunn et. al. 2008).
Summation is over j=1 to N, j ≠ a. Total is N-1.

Overall mean mutual information (Dunn et. al. 2008).
2/(N*(N-1)) by the sum of MI where the indices run i=1 to N-1, j=i+1 to N (triu).
"""
function _mean_total(mi::Matrix)
  values = matrix2list(mi)
  sum(values) / length(values)
end

"""
APC (Dunn et. al. 2008)
"""
function APC!{T}(MI::Matrix{T})
  nrow, ncol = size(MI)
  MI_mean = _mean_total(MI)
#   if MI_mean ==  0.0
#     return(fill!(MI, zero(T)))
#   end
  MI_res_mean = _mean_column(MI)
  @inbounds for j in 1:ncol
    for i in 1:nrow
      if i != j
        MI[i,j] -= ((MI_res_mean[i] * MI_res_mean[j]) /  MI_mean)
      else
        MI[i,j] = NaN
      end
    end
  end
  MI
end
