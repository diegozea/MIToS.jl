# Pseudofrequencies
# =================
#
# BLOSUM based pseudofrequencies
#

"""
`blosum_pseudofrequencies!{T}(Gab::ContingencyTable{T,2,UngappedAlphabet}, Pab::ContingencyTable{T,2,UngappedAlphabet})`

This function uses the conditional probability matrix `BLOSUM62_Pij` to fill a preallocated
`Gab` with pseudo frequencies. `blosum_pseudofrequencies!` also needs the real
frequencies/probabilities `Pab`. This observed probabilities are then used to estimate
the pseudo frequencies.

`Gab = Σcd  Pcd ⋅ BLOSUM62( a | c ) ⋅ BLOSUM62( b | d )`
"""
function blosum_pseudofrequencies!{T}(Gab::ContingencyTable{T,2,UngappedAlphabet},
                                      Pab::ContingencyTable{T,2,UngappedAlphabet})
    @assert get_total(Pab) ≈ one(T) "The input should be a probability table (normalized)"
    pab  = array(get_table(Pab))
    gab  = array(get_table(Gab))
    bl62 = array(get_table(BLOSUM62_Pij))
    @inbounds for a in 1:20, b in 1:20
        gab[a,b] = zero(T)
        for i in 1:20
            bl_ia = bl62[i,a]
            for j in 1:20
                P = pab[i,j]
                if P != 0
                    # BLOSUM62_P_i_j[i,a] is p(a | i)
                    gab[a,b] += ( P * bl_ia * bl62[j,b] )
                end
            end
        end
    end
    update_marginals!(Gab)
    normalize!(Gab)
end

"""
`apply_pseudofrequencies!{T}(Pab::ContingencyTable{T,2,UngappedAlphabet}, Gab::ContingencyTable{T,2,UngappedAlphabet}, α, β)`

Apply pseudofrequencies `Gab` over `Pab`, as a weighted mean of both. α is the weight of
the real frequencies `Pab` and β the weight of the pseudofrequencies.

`Pab = (α ⋅ Pab + β ⋅ Gab )/(α + β)`
"""
function apply_pseudofrequencies!{T}(Pab::ContingencyTable{T,2,UngappedAlphabet},
                                     Gab::ContingencyTable{T,2,UngappedAlphabet}, α, β)
    @assert get_total(Pab) ≈ one(T) "Pab should be a probability table (normalized)"
    @assert get_total(Gab) ≈ one(T) "Gab should be a probability table (normalized)"
    if β == 0.0
        return(Pab)
    end
    pab  = array(get_table(Pab))
    gab  = array(get_table(Gab))
    frac = one(T) / ( α + β )
    @inbounds for i in 1:20
        @simd for j in 1:20
            pab[i,j] = ( α * pab[i,j] + β * gab[i,j] ) * frac
        end
    end
    update_marginals!(Pab)
    normalize!(Pab)
end
