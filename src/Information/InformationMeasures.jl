# Entropy
# =======

"""
    shannon_entropy(table::Union{Counts{T,N,A},Probabilities{T,N,A}}; base::Union{Real, Nothing}=nothing)

It calculates the Shannon entropy (H) from a table of `Counts` or `Probabilities`.
Use last and optional positional argument to change the base of the log. The default base
(`base=nothing`) is ℯ, so the result is in nats. You can use `base=2.0` to get the result 
in bits.
"""
function shannon_entropy(table::Probabilities{T,N,A}; base::Union{Real, Nothing}=nothing) where {T,N,A}
    H = zero(T)
    p = gettablearray(table)
    @inbounds for pᵢ in p
        if pᵢ > zero(T)
            H -= pᵢ * log(pᵢ)
        end
    end
    base === nothing ? H : (H / log(base))
end

function shannon_entropy(table::Counts{T,N,A}; base::Union{Real, Nothing}=nothing) where {T,N,A}
    H = zero(T)
    total = gettotal(table)
    n = gettablearray(table)
    @inbounds for nᵢ in n
        if nᵢ > zero(T)
            H -= nᵢ * log(nᵢ / total)
        end
    end
    if base === nothing
        H / total  # Default base: e
    else
        (H / total) / log(base)
    end
end

function StatsBase.entropy(table::Union{Counts{T,N,A},Probabilities{T,N,A}}) where {T,N,A}
    Base.depwarn("entropy(table::Union{Counts,Probabilities}) is deprecated. Use shannon_entropy(table) instead.", :entropy, force=true)
    shannon_entropy(table)
end

function StatsBase.entropy(table::Union{Counts{T,N,A},Probabilities{T,N,A}}, base::Real) where {T,N,A}
    Base.depwarn("entropy(table::Union{Counts,Probabilities}, base::Real) is deprecated. Use shannon_entropy(table; base=base) instead.", :entropy, force=true)
    shannon_entropy(table, base=base)
end

# using mapfreq to define the method for multiple sequence alignments
"""
    shannon_entropy(msa::AbstractArray{Residue}; base::Union{Real, Nothing}=nothing, 
        probabilities::Bool=false, usediagonal::Bool=true, kargs...)

It calculates the Shannon entropy (H) on a MSA. You can use the keyword argument `base` to
change the base of the log. The default base is ℯ (`base=nothing`), so the result is in nats. 
You can use 2.0 as base to get the result in bits. It uses [`mapfreq`](@ref) under the hood, 
so it takes the same keyword arguments. By default, it measures the entropy of each column 
in the MSA. You can use `dims = 1` to measure the entropy of each sequence. You can also 
set `rank = 2`to measure the joint entropy of each pair of sequences or columns. This 
function sets by default the `probabilities` keyword argument to `false` because it's 
faster to calculate the entropy from counts/frequencies. It also sets `usediagonal = true`
to also calculate the entropy of the individual variables (sequences or columns).

```jldoctest
julia> using MIToS.MSA, MIToS.Information

julia> msa = Residue['C' 'G'; 'C' 'L'; 'C' 'I']
3×2 Matrix{Residue}:
 C  G
 C  L
 C  I

julia> shannon_entropy(msa)
1×2 Named Matrix{Float64}
 Function ╲ Col │       1        2
────────────────┼─────────────────
shannon_entropy │     0.0  1.09861

"""
function shannon_entropy(msa::AbstractArray{Residue}; probabilities::Bool=false, 
    usediagonal=true, kargs...)
    mapfreq(shannon_entropy, msa; probabilities=probabilities, usediagonal=usediagonal, kargs...)
end

# Marginal Entropy
# ----------------

function _marginal_entropy(table::Probabilities{T,N,A}, margin::Int) where {T,N,A}
    H = zero(T)
    marginals = getmarginalsarray(table)
    @inbounds for pi in view(marginals, :, margin)
        if pi > zero(T)
            H += pi * log(pi)
        end
    end
    -H # Default base: e
end

function _marginal_entropy(table::Counts{T,N,A}, margin::Int) where {T,N,A}
    H = zero(T)
    total = gettotal(table)
    marginals = getmarginalsarray(table)
    @inbounds for ni in view(marginals, :, margin)
        if ni > zero(T)
            H += ni * log(ni / total)
        end
    end
    -H / total # Default base: e
end

"""
    marginal_entropy(table::Union{Counts{T,N,A},Probabilities{T,N,A}}; margin::Int=1, 
        base::Union{Real, Nothing}=nothing)

It calculates marginal entropy (H) from a table of `Counts` or `Probabilities`. It takes 
two keyword arguments: `margin` and `base`. The first one is used to indicate the margin
used to calculate the entropy, e.g. it estimates the entropy H(X) if margin is 1, H(Y)
for 2, etc. The default value of `margin` is 1. The second keyword argument is used to
change the base of the log. The default base is ℯ (`base = nothing`), so the result is in
nats. You can use `base = 2.0` to get the result in bits.
"""
function marginal_entropy(
    table::Union{Counts{T,N,A},Probabilities{T,N,A}};
    margin::Int=1,
    base::Union{Real, Nothing}=nothing,
) where {T,N,A}
    H = _marginal_entropy(table, margin)
    if base === nothing
        H # Default base: e
    else
        H / log(base)
    end
end

# Deprecate the marginal_entropy methods taking positional arguments
function marginal_entropy(table::Union{Counts{T,N,A},Probabilities{T,N,A}}, margin::Int) where {T,N,A}
    Base.depwarn("marginal_entropy(table, margin) is deprecated. Use marginal_entropy(table; margin=margin) instead.", :marginal_entropy, force=true)
    marginal_entropy(table, margin=margin)
end

function marginal_entropy(table::Union{Counts{T,N,A},Probabilities{T,N,A}}, margin::Int, base::Real) where {T,N,A}
    Base.depwarn("marginal_entropy(table, margin, base) is deprecated. Use marginal_entropy(table; margin=margin, base=base) instead.", :marginal_entropy, force=true)
    marginal_entropy(table, margin=margin, base=base)
end

# Kullback-Leibler
# ================

function _gettablearray(table::Union{Probabilities{T,N,A}, Counts{T,N,A}, ContingencyTable{T,N,A}}) where {T,N,A}
    gettablearray(table)
end
_gettablearray(table::Array{T,N}) where {T,N} = table

const KL_KARG_DOC = """
You can use the keyword argument `background` to set the background distribution. This 
argument can take an `Array`, `Probabilities`, or `ContingencyTable` object. The background 
distribution must have the same size and alphabet as the probabilities. The default is the 
`BLOSUM62_Pi` table.  You can use the `base` keyword argument to change the base of the log.
The default base of the log is ℯ (`base = nothing`).
"""

"""
    kullback_leibler(probabilities::Probabilities{T,N,A}, background::Union{Array{T,N}, 
        Probabilities{T,N,A}, ContingencyTable{T,N,A}}=BLOSUM62_Pi, 
        base::Union{Real, Nothing}=nothing)

It calculates the Kullback-Leibler (KL) divergence from a table of `Probabilities`. 
$KL_KARG_DOC
"""
function kullback_leibler(
    probabilities::Probabilities{T,N,A};
    background::Union{Array{T,N}, Probabilities{T,N,A}, ContingencyTable{T,N,A}} = BLOSUM62_Pi,
    base::Union{Real, Nothing}=nothing) where {T,N,A}
    p = getcontingencytable(probabilities)
    bg = _gettablearray(background)
    @assert size(background) == size(p) "probabilities and background must have the same size."
    KL = zero(T)
    @inbounds for i in eachindex(p)
        pᵢ = p[i]
        if pᵢ > zero(T)
            KL += pᵢ * log(pᵢ / bg[i])
        end
    end
    if base === nothing
        KL # Default base: e
    else
        KL / log(base)
    end
end

# Kullback-Leibler for MSA

"""
    kullback_leibler(msa::AbstractArray{Residue}; background::Union{Array{T,N}, Probabilities{T,N,A}, ContingencyTable{T,N,A}}=BLOSUM62_Pi, base::Union{Real, Nothing}=nothing, kargs...)

It calculates the Kullback-Leibler (KL) divergence from a multiple sequence alignment (MSA).
$KL_KARG_DOC The other keyword arguments are passed to the [`mapfreq`](@ref) function.
"""
function kullback_leibler(msa::AbstractArray{Residue}; background::Union{Array{T,N}, Probabilities{T,N,A}, ContingencyTable{T,N,A}}=BLOSUM62_Pi, base::Union{Real, Nothing}=nothing, rank::Int=1, kargs...) where {T,N,A}
    @assert rank == 1 "rank must be 1 for kullback_leibler"
    mapfreq(kullback_leibler, msa; rank=rank, background=background, base=base, kargs...)
end

# Deprecate the old methods

# Method with positional arguments for background and base
function kullback_leibler(p::Probabilities{T,N,A}, q::Union{Array{T,N}, Probabilities{T,N,A}, ContingencyTable{T,N,A}}, base::Real) where {T,N,A}
    Base.depwarn("kullback_leibler(p, q, base) is deprecated. Use kullback_leibler(p; background=q, base=base) instead.", :kullback_leibler, force=true)
    kullback_leibler(p; background=q, base=base)
end

# Method with positional argument for background
function kullback_leibler(p::Probabilities{T,N,A}, q::Union{Array{T,N}, Probabilities{T,N,A}, ContingencyTable{T,N,A}}) where {T,N,A}
    Base.depwarn("kullback_leibler(p, q) is deprecated. Use kullback_leibler(p; background=q) instead.", :kullback_leibler, force=true)
    kullback_leibler(p; background=q)
end

# Method with positional argument for base
function kullback_leibler(p::Probabilities{T,N,A}, base::Real) where {T,N,A}
    Base.depwarn("kullback_leibler(p, base) is deprecated. Use kullback_leibler(p; base=base) instead.", :kullback_leibler, force=true)
    kullback_leibler(p; base=base)
end

# Mutual Information
# ==================

# It avoids ifelse() because log is expensive (https://github.com/JuliaLang/julia/issues/8869)
@inline function _mi(::Type{T}, pij, pi, pj) where {T}
    (pij > zero(T)) && (pi > zero(T)) ? T(pij * log(pij / (pi * pj))) : zero(T)
end

function mutual_information(table::Probabilities{T,2,A}) where {T,A}
    MI = zero(T)
    marginals = getmarginalsarray(table)
    p = gettablearray(table)
    N = size(marginals, 1)
    @inbounds for j = 1:N
        pj = marginals[j, 2]
        if pj > zero(T)
            @inbounds @simd for i = 1:N
                MI += _mi(T, p[i, j], marginals[i, 1], pj)
            end
        end
    end
    MI # Default base: e
end

# It avoids ifelse() because log is expensive (https://github.com/JuliaLang/julia/issues/8869)
@inline function _mi(total::T, nij, ni, nj) where {T}
    (nij > zero(T)) && (ni > zero(T)) ? T(nij * log((total * nij) / (ni * nj))) : zero(T)
end

function mutual_information(table::Counts{T,2,A}) where {T,A}
    MI = zero(T)
    marginals = getmarginalsarray(table)
    n = gettablearray(table)
    total = gettotal(table)
    N = size(marginals, 1)
    @inbounds for j = 1:N
        nj = marginals[j, 2]
        if nj > zero(T)
            @inbounds @simd for i = 1:N
                MI += _mi(total, n[i, j], marginals[i, 1], nj)
            end
        end
    end
    MI / total # Default base: e
end

"""
It calculates Mutual Information (MI) from a table of `Counts` or `Probabilities`.
Use last and optional positional argument to change the base of the log. The default base
is e, so the result is in nats. You can use 2.0 as base to get the result in bits.
Calculation of MI from `Counts` is faster than from `Probabilities`.
"""
function mutual_information(
    table::Union{Counts{T,N,A},Probabilities{T,N,A}},
    base::Real,
) where {T,N,A}
    mutual_information(table) / log(base)
end

function mutual_information(pxyz::Union{Counts{T,3,A},Probabilities{T,3,A}}) where {T,A}
    pxy = delete_dimensions(pxyz, 3)
    return (
        marginal_entropy(pxyz, margin=1) +                 # H(X) +
        marginal_entropy(pxyz, margin=2) +                 # H(Y) +
        marginal_entropy(pxyz, margin=3) -                 # H(Z) -
        shannon_entropy(pxy) -                              # H(X,Y) -
        shannon_entropy(delete_dimensions!(pxy, pxyz, 2)) - # H(X,Z) -
        shannon_entropy(delete_dimensions!(pxy, pxyz, 2)) + # H(Y,Z) +
        shannon_entropy(pxyz)                               # H(X,Y,Z)
    )
end

# Normalized Mutual Information by Entropy
# ----------------------------------------

"""
It calculates a Normalized Mutual Information (nMI) by Entropy from a table of `Counts` or
`Probabilities`.

`nMI(X, Y) = MI(X, Y) / H(X, Y)`
"""
function normalized_mutual_information(
    table::Union{Counts{T,N,A},Probabilities{T,N,A}},
) where {T,N,A}
    H = shannon_entropy(table)
    if H != zero(T)
        MI = mutual_information(table)
        return (T(MI / H))
    else
        return (zero(T))
    end
end
