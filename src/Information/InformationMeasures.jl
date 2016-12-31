# Entropy
# =======

# """
# Shannon entropy (H)
# """
# """
# `estimate(Entropy(base), p)`
#
# `p` should be a `ResidueProbability` table. The result type is determined by `base`.
# """
function entropy{T,N,A}(table::Probabilities{T,N,A})
    H = zero(T)
    p = gettablearray(table)
    @inbounds for pᵢ in p
        if pᵢ != zero(T)
            H -= pᵢ * log(pᵢ)
        end
    end
    H # Default base: e
end

# """
# `estimate(Entropy{T}(base), n::ResidueCount)`
#
# It's the fastest option (you don't spend time on probability calculations).
# The result type is determined by the `base`.
# """
function entropy{T,N,A}(table::Counts{T,N,A})
    H = zero(T)
    total = gettotal(table)
    n = gettablearray(table)
    @inbounds for nᵢ in n
        if nᵢ != zero(T)
            H -= nᵢ * log(nᵢ/total)
        end
    end
    H/total # Default base: e
end

function entropy{T,N,A}(table::Union{Counts{T,N,A},Probabilities{T,N,A}}, base::T)
    entropy(table) / log(base)
end

# Marginal Entropy
# ----------------

# """
# `estimate_on_marginal(Entropy{T}(base), p, marginal)`
#
# This function estimate the entropy H(X) if marginal is 1, H(Y) for 2, etc.
# The result type is determined by `base`.
# """
function marginal_entropy{T,N,A}(table::Probabilities{T,N,A}, margin::Int)
    H = zero(T)
    marginals = getmarginalsarray(table)
    @inbounds for pi in view(marginals, :, margin)
        if pi != zero(T)
            H -= pi * log(pi)
        end
    end
    H # Default base: e
end

function marginal_entropy{T,N,A}(table::Counts{T,N,A}, margin::Int)
    H = zero(T)
    total = gettotal(table)
    marginals = getmarginalsarray(table)
    @inbounds for ni in view(marginals, :, margin)
        if ni != zero(T)
            H -= ni * log(ni/total)
        end
    end
    H/total # Default base: e
end

function marginal_entropy{T,N,A}(table::Union{Counts{T,N,A},Probabilities{T,N,A}},
                                 margin::Int, base::T)
    marginal_entropy(table, margin) / log(base)
end

# Kullback-Leibler
# ================


# """
# Kullback-Leibler (KL). This `SymmetricMeasure` has two fields.
# The first is the base of the logarithm and the second is the backgroud frequency.
# """
# """
# `estimate(KullbackLeibler(base, background), p)`
#
# `p` should be a `ResidueProbability` table, and `background` must have the size of `p`.
# The result type is determined by `base` and `background`.
# """
function kullback_leibler{T,A}(probabilities::Probabilities{T,1,A},
                               background::Vector{T})
    p = getcontingencytable(probabilities)
    @assert size(background)!=size(p) "probabilities and background must have the same size."
    KL = zero(T)
    @inbounds for i in 1:length(p)
        pi = p[i]
        if pi != 0.0
            KL += pi * log(pi/background[i])
        end
    end
    KL # Default base: e
end

function kullback_leibler{T,A}(probabilities::Probabilities{T,1,A},
                               background::ContingencyTable{T,1,A}=BLOSUM62_Pi)
    kullback_leibler(probabilities, gettablearray(background))
end

kullback_leibler{T,A}(p::Probabilities{T,1,A}, q, base::T) = kullback_leibler(p, q)/log(base)
kullback_leibler{T,A}(p::Probabilities{T,1,A}, base::T) = kullback_leibler(p)/log(base)

# Mutual Information
# ==================

# It avoids ifelse() because log is expensive (https://github.com/JuliaLang/julia/issues/8869)
@inline function _mi{T}(::Type{T}, pij, pi, pj)
    @fastmath (pij > zero(T)) && (pi > zero(T)) ? T(pij * log(pij/(pi*pj))) : zero(T)
end

# """
# Mutual Information (MI)
# """
# """
# `estimate(MutualInformation(), pxy::ResidueProbability [, base])`
#
# Calculate Mutual Information from `ResidueProbability`. The result type is determined by `base`.
# """
function mutual_information{T,A}(table::Probabilities{T,2,A})
    MI = zero(T)
    marginals = getmarginalsarray(table)
    p = gettablearray(table)
    N = size(marginals,1)
    @inbounds for j in 1:N
        pj = marginals[j,2]
        if pj > 0.0
            @inbounds @simd for i in 1:N
                MI += _mi(T, p[i,j], marginals[i,1], pj)
            end
        end
    end
    MI # Default base: e
end

# It avoids ifelse() because log is expensive (https://github.com/JuliaLang/julia/issues/8869)
@inline function _mi{T}(total::T, nij, ni, nj)
    @fastmath (nij > zero(T)) && (ni > zero(T)) ? T(nij * log((total * nij)/(ni * nj))) : zero(T)
end

# """
# `estimate(MutualInformation(), pxy::ResidueCount [, base])`
#
# Calculate Mutual Information from `ResidueCount`. The result type is determined by the `base`.
# It's the fastest option (you don't spend time on probability calculations).
# """
function mutual_information{T,A}(table::Counts{T,2,A})
    MI = zero(T)
    marginals = getmarginalsarray(table)
    n = gettablearray(table)
    total = gettotal(table)
    N = size(marginals,1)
    @inbounds for j in 1:N
        nj = marginals[j,2]
        if nj > 0.0
            @inbounds @simd for i in 1:N
                MI += _mi(total, n[i,j], marginals[i,1], nj)
            end
        end
    end
    MI/total # Default base: e
end

function mutual_information{T,N,A}(table::Union{Counts{T,N,A}, Probabilities{T,N,A}},
                                   base::T)
    mutual_information(table) / log(base)
end

function mutual_information{T,A}(pxyz::Probabilities{T,3,A})
    pxy = delete_dimensions(pxyz, 3)
    return(
        marginal_entropy(pxyz, 1) +                 # H(X) +
        marginal_entropy(pxyz, 2) +                 # H(Y) +
        marginal_entropy(pxyz, 3) -                 # H(Z) -
        entropy(pxy) -                              # H(X,Y) -
        entropy(delete_dimensions!(pxy, pxyz, 2)) - # H(X,Z) -
        entropy(delete_dimensions!(pxy, pxyz, 2)) + # H(Y,Z) +
        entropy(pxyz)                               # H(X,Y,Z)
    )
end

# Normalized Mutual Information by Entropy
# ----------------------------------------

# """
# Normalized Mutual Information (nMI) by Entropy.
#
# `nMI(X, Y) = MI(X, Y) / H(X, Y)`
# """
function normalized_mutual_information{T,N,A}(table::Union{Counts{T,N,A},Probabilities{T,N,A}})
    H = entropy(table)
    if H != zero(T)
        MI = mutual_information(table)
        return(T(MI/H))
    else
        return(zero(T))
    end
end

# Pairwise Gap Percentage
# =======================

function gap_intersection_percentage{T}(nxy::Counts{T,2,GappedAlphabet})
    T(100.0) * nxy[21,21] / gettotal(nxy)
end

function gap_union_percentage{T}(nxy::Counts{T,2,GappedAlphabet})
    marginals = getmarginalsarray(nxy)
    T(100.0) * (marginals[21,1] + marginals[21,2] - nxy[21,21]) / gettotal(nxy)
end
