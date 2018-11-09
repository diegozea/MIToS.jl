"""
Mean mutual information of column a (Dunn et. al. 2008).
Summation is over j=1 to N, j ≠ a. Total is N-1.
"""
function _mean_column(mi::Matrix{T}) where {T}
    (dropdims(sum(mi, dims=1), dims=1) .- diag(mi)) ./ (size(mi, 1)-1)
end

"""
Mean mutual information of column a (Dunn et. al. 2008).
Summation is over j=1 to N, j ≠ a. Total is N-1.

Overall mean mutual information (Dunn et. al. 2008).
2/(N*(N-1)) by the sum of MI where the indices run i=1 to N-1, j=i+1 to N (triu).
"""
function _mean_total(mi::Matrix{T}) where T
    values = matrix2list(mi)
    sum(values) / length(values)
end

"APC (Dunn et. al. 2008)"
function APC!(MI::Matrix{T}) where T
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

function APC!(MI::PairwiseListMatrix{T,true,Vector{T}}) where T <: AbstractFloat
    ncol = MI.nelements
    MI_col_sum = vec(sum_nodiag(MI, dims=1))
    MI_mean = sum(MI_col_sum) / (length(MI) - ncol)
    MI_col_mean = MI_col_sum ./ (ncol-1)
    k = 0
    list = MI.list
    @inbounds for i in 1:ncol
        for j in i:ncol
            k += 1
            if i != j
                list[k] -= ((MI_col_mean[i] * MI_col_mean[j]) /  MI_mean)
            else
                list[k] = NaN
            end
        end
    end
    MI
end

function APC!(MI::PairwiseListMatrix{T,false,Vector{T}}) where T <: AbstractFloat
    ncol = MI.nelements
    MI_col_sum = vec(sum_nodiag(MI, dims=1))
    MI_mean = sum(MI_col_sum) / (length(MI) - ncol)
    MI_col_mean = MI_col_sum ./ (ncol-1)
    k = 0
    list = MI.list
    @inbounds for i in 1:(ncol-1)
        for j in (i+1):ncol
            k += 1
            list[k] -= ((MI_col_mean[i] * MI_col_mean[j]) /  MI_mean)
        end
    end
    MI.diag[:] .= NaN
    MI
end

function APC!(MI::NamedArray{T,2,PairwiseListMatrix{T,D,Vector{T}},NTuple{2,OrderedDict{String,Int}}}) where {T,D}
    plm = getarray(MI)::PairwiseListMatrix{T,D,Vector{T}}
    APC!(plm)
    MI
end
