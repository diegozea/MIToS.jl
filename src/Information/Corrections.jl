"""
APC (Dunn et. al. 2008)
"""
function APC!(MI::Matrix)
  nrow, ncol = size(MI)
  MI_mean = mean(MI)
  MI_res_mean = mean(MI, 1)
  @inbounds for j in 1:ncol
    for i in 1:nrow
      MI[i,j] -= ((MI_res_mean[i] * MI_res_mean[j]) /  MI_mean)
    end
  end
  MI
end
