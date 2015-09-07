"""
APC (Dunn et. al. 2008)
"""
function APC!{T}(MI::Matrix{T}; diagonal::Bool=true)
  nrow, ncol = size(MI)
  MI_mean = mean(MI)
  if MI_mean ==  0.0
    return(fill!(MI, zero(T)))
  end
  MI_res_mean = mean(MI, 1)
  @inbounds for j in 1:ncol
    for i in 1:nrow
      if i != j || diagonal
        MI[i,j] -= ((MI_res_mean[i] * MI_res_mean[j]) /  MI_mean)
      end
    end
  end
  MI
end
