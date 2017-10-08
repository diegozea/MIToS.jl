# Entropy
# =======

function StatsBase.entropy(table::Probabilities{T,N,A}) where {T,N,A}
    H = zero(T)
    p = gettablearray(table)
    @inbounds for pᵢ in p
        if pᵢ > zero(T)
            H += pᵢ * log(pᵢ)
        end
    end
    -H # Default base: e
end

function StatsBase.entropy(table::Counts{T,N,A}) where {T,N,A}
    H = zero(T)
    total = gettotal(table)
    n = gettablearray(table)
    @inbounds for nᵢ in n
        if nᵢ > zero(T)
            H += nᵢ * log(nᵢ/total)
        end
    end
    -H/total # Default base: e
end

"""
It calculates the Shannon entropy (H) from a table of `Counts` or `Probabilities`.
Use last and optional positional argument to change the base of the log. The default base
is e, so the result is in nats. You can use 2.0 as base to get the result in bits.
"""
function StatsBase.entropy(table::Union{Counts{T,N,A},Probabilities{T,N,A}}, base::Real) where {T,N,A}
    entropy(table) / log(base)
end

# Marginal Entropy
# ----------------

function marginal_entropy(table::Probabilities{T,N,A}, margin::Int) where {T,N,A}
    H = zero(T)
    marginals = getmarginalsarray(table)
    @inbounds for pi in view(marginals, :, margin)
        if pi > zero(T)
            H += pi * log(pi)
        end
    end
    -H # Default base: e
end

function marginal_entropy(table::Counts{T,N,A}, margin::Int) where {T,N,A}
    H = zero(T)
    total = gettotal(table)
    marginals = getmarginalsarray(table)
    @inbounds for ni in view(marginals, :, margin)
        if ni > zero(T)
            H += ni * log(ni/total)
        end
    end
    -H/total # Default base: e
end

"""
It calculates marginal entropy (H) from a table of `Counts` or `Probabilities`. The second
positional argument is used to indicate the magin used to calculate the entropy, e.g. it
estimates the entropy H(X) if marginal is 1, H(Y) for 2, etc.
Use last and optional positional argument to change the base of the log. The default base
is e, so the result is in nats. You can use 2.0 as base to get the result in bits.
"""
function marginal_entropy(table::Union{Counts{T,N,A},Probabilities{T,N,A}},
                          margin::Int, base::Real) where {T,N,A}
    marginal_entropy(table, margin) / log(base)
end

# Kullback-Leibler
# ================


function kullback_leibler(probabilities::Probabilities{T,N,A},
                          background::Array{T,N}) where {T,N,A}
    p = getcontingencytable(probabilities)
    @assert size(background)==size(p) "probabilities and background must have the same size."
    KL = zero(T)
    @inbounds for i in 1:length(p)
        pi = p[i]
        if pi > zero(T)
            KL += pi * log(pi/background[i])
        end
    end
    KL # Default base: e
end

function kullback_leibler(probabilities::Probabilities{T,N,A},
         background::Union{Probabilities{T,N,A},ContingencyTable{T,N,A}}=BLOSUM62_Pi) where {T,N,A}
    kullback_leibler(probabilities, gettablearray(background))
end

"""
It calculates the Kullback-Leibler (KL) divergence from a table of `Probabilities`. The
second positional argument is a `Probabilities` or `ContingencyTable` with the background
distribution. It's optional, the default is the `BLOSUM62_Pi` table.
Use last and optional positional argument to change the base of the log. The default base
is e, so the result is in nats. You can use 2.0 as base to get the result in bits.
"""
kullback_leibler(p::Probabilities{T,1,A}, q, base::Real) where {T,A} = kullback_leibler(p, q)/log(base)
kullback_leibler(p::Probabilities{T,1,A}, base::Real) where {T,A} = kullback_leibler(p)/log(base)

# Mutual Information
# ==================

# It avoids ifelse() because log is expensive (https://github.com/JuliaLang/julia/issues/8869)
@inline function _mi(::Type{T}, pij, pi, pj) where T
    (pij > zero(T)) && (pi > zero(T)) ? T(pij * log(pij/(pi*pj))) : zero(T)
end

function mutual_information(table::Probabilities{T,2,A}) where {T,A}
    MI = zero(T)
    marginals = getmarginalsarray(table)
    p = gettablearray(table)
    N = size(marginals,1)
    @inbounds for j in 1:N
        pj = marginals[j,2]
        if pj > zero(T)
            @inbounds @simd for i in 1:N
                MI += _mi(T, p[i,j], marginals[i,1], pj)
            end
        end
    end
    MI # Default base: e
end

# It avoids ifelse() because log is expensive (https://github.com/JuliaLang/julia/issues/8869)
@inline function _mi(total::T, nij, ni, nj) where T
    (nij > zero(T)) && (ni > zero(T)) ? T(nij * log((total * nij)/(ni * nj))) : zero(T)
end

function mutual_information(table::Counts{T,2,A}) where {T,A}
    MI = zero(T)
    marginals = getmarginalsarray(table)
    n = gettablearray(table)
    total = gettotal(table)
    N = size(marginals,1)
    @inbounds for j in 1:N
        nj = marginals[j,2]
        if nj > zero(T)
            @inbounds @simd for i in 1:N
                MI += _mi(total, n[i,j], marginals[i,1], nj)
            end
        end
    end
    MI/total # Default base: e
end

"""
It calculates Mutual Information (MI) from a table of `Counts` or `Probabilities`.
Use last and optional positional argument to change the base of the log. The default base
is e, so the result is in nats. You can use 2.0 as base to get the result in bits.
Calculation of MI from `Counts` is faster than from `Probabilities`.
"""
function mutual_information(table::Union{Counts{T,N,A}, Probabilities{T,N,A}},
                            base::Real) where {T,N,A}
    mutual_information(table) / log(base)
end

function mutual_information(pxyz::Union{Counts{T,3,A}, Probabilities{T,3,A}}) where {T,A}
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

"""
It calculates a Normalized Mutual Information (nMI) by Entropy from a table of `Counts` or
`Probabilities`.

`nMI(X, Y) = MI(X, Y) / H(X, Y)`
"""
function normalized_mutual_information(table::Union{Counts{T,N,A},Probabilities{T,N,A}}) where {T,N,A}
    H = entropy(table)
    if H != zero(T)
        MI = mutual_information(table)
        return(T(MI/H))
    else
        return(zero(T))
    end
end
