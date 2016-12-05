# Pseudofrequencies
# =================

abstract Pseudofrequencies

immutable NoPseudofrequencies <: Pseudofrequencies end

# BLOSUM based pseudofrequencies
# ==============================

immutable BLOSUM_Pseudofrequencies
    α::Float64
    β::Float64
end

"""
`_calculate_blosum_pseudofrequencies!{T}(Pab::ContingencyTable{T,2,UngappedAlphabet})`

This function uses the conditional probability matrix `BLOSUM62_Pij` to fill the temporal
array field of `Pab` with pseudo frequencies (`Gab`). This function needs the real
frequencies/probabilities `Pab` because they are used to estimate the pseudofrequencies.

`Gab = Σcd  Pcd ⋅ BLOSUM62( a | c ) ⋅ BLOSUM62( b | d )`
"""
function _calculate_blosum_pseudofrequencies!{T}(Pab::ContingencyTable{T,2,UngappedAlphabet})
    @assert get_total(Pab) ≈ one(T) "The input should be a probability table (normalized)"
    pab  = array(get_table(Pab))
    gab  = Pab.temporal
    bl62 = array(get_table(BLOSUM62_Pij))
    total = zero(T)
    @inbounds for b in 1:20, a in 1:20
        gab[a,b] = zero(T)
        for d in 1:20
            bl62_db = bl62[d,b]
            for c in 1:20
                P = pab[c,d]
                if P != 0
                    # BLOSUM62_Pij[c,a] is p(a|c)
                    gab[a,b] += ( P * bl62[c,a] * bl62_db )
                end
            end
        end
        total += g[a,b]
    end
    if total ≉ one(T)
        @inbounds for col in 1:20
            @simd for row in 1:20
                gab[row,col] /= total
            end
        end
    end
    Pab
end

"""
`apply_pseudofrequencies!{T}(Pab::ContingencyTable{T,2,UngappedAlphabet}, pseudofrequencies::BLOSUM_Pseudofrequencies)`

When a `BLOSUM_Pseudofrequencies(α,β)` is used, this function applies pseudofrequencies
`Gab` over `Pab`, as a weighted mean of both. It uses the conditional probability
matrix `BLOSUM62_Pij` and the real frequencies/probabilities `Pab` to estimate the
pseudofrequencies `Gab`. α is the weight of the real frequencies `Pab` and β the weight
of the pseudofrequencies.

`Gab = Σcd  Pcd ⋅ BLOSUM62( a | c ) ⋅ BLOSUM62( b | d )`
`Pab = (α ⋅ Pab + β ⋅ Gab )/(α + β)`
"""
function apply_pseudofrequencies!{T}(Pab::ContingencyTable{T,2,UngappedAlphabet},
                                     pseudofrequencies::BLOSUM_Pseudofrequencies)
    α = T(pseudofrequencies.α)
    β = T(pseudofrequencies.β)
    if β == 0.0
        return(Pab)
    end
    _calculate_blosum_pseudofrequencies!(Pab)
    pab  = array(get_table(Pab))
    gab  = Pab.temporal
    frac = one(T) / ( α + β )
    @inbounds for col in 1:20
        @simd for row in 1:20
            pab[row,col] = ( α * pab[row,col] + β * gab[row,col] ) * frac
        end
    end
    update_marginals!(Pab)
    normalize!(Pab)
end

@inline apply_pseudofrequencies!(Pab::ContingencyTable, pseudofrequencies::NoPseudofrequencies) = Pab
