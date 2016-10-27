# Calculates percent identity of two aligned sequences
# No account of the positions with gaps in both sequences in the length

"seq1 and seq2 should have the same len"
function _percentidentity(seq1, seq2, len)
    count = 0
    colgap = 0
    colxaa = 0
    @inbounds for i in 1:len
        aa1 = seq1[i]
        aa2 = seq2[i]
        # Columns with Residue('X') aren't used
        if aa1 == XAA || aa2 == XAA
            colxaa += 1
            continue
        end
        if aa1 == aa2
            count += 1
            colgap += Int(aa1 == GAP)
        end
    end
    100.0 * (count-colgap)/(len-colgap-colxaa)
end

"""
`percentidentity(seq1, seq2)`

Calculates the fraction of identities between two aligned sequences. The identity value is
calculated as the number of identical characters in the i-th position of both sequences
divided by the length of both sequences. Positions with gaps in both sequences doesn't
count to the length of the sequences. Positions with a `XAA` in at least one sequence
aren't counted.
"""
function percentidentity(seq1, seq2)
    len = length(seq1)
    if len != length(seq2)
        throw(ErrorException("""
        Sequences of different length, they aren't aligned or don't come from the same MSA.
        """))
    end
    _percentidentity(seq1, seq2, len)
end

"""
`percentidentity(seq1, seq2, threshold)`

Computes quickly if two aligned sequences have a identity value greater than a given
`threshold` value. Returns a boolean value. Positions with gaps in both sequences
doesn't count to the length of the sequences. Positions with a `XAA` in at least one
sequence aren't counted.
"""
function percentidentity(seq1, seq2, threshold)
    fraction = threshold / 100.0
    len = length(seq1)
    if len != length(seq2)
        throw(ErrorException("""
        Sequences of different length, they aren't aligned or don't come from the same MSA.
        """))
    end
    n = len
    limit_count = n * fraction
    diff = 0
    count = 0
    @inbounds for i in 1:len
        aa1 = seq1[i]
        aa2 = seq2[i]
        if aa1 == XAA || aa2 == XAA
            n -= 1
            continue
        end
        if aa1 == aa1
            if aa1 != GAP
                count += 1
                if count >= limit_count
                    return(true)
                end
            else
                n -= 1
                limit_count = n * fraction
            end
        else
            diff += 1
            if diff > n - limit_count
                return(false)
            end
        end
    end
    (count/n) >= fraction
end

# percentidentity for a MSA
# =========================

function _percentidentity_kernel!(scores, aln::Vector{Vector{Residue}}, nseq, len)
    k = 0
    list = getlist(scores)
    @inbounds for i in 1:(nseq-1)
        a = aln[i]
        for j in (i+1):nseq
            list[k += 1] = _percentidentity(a, aln[j], len)
        end
    end
    scores
end

"""
`percentidentity(msa[, out::Type=Float64])`

Calculates the identity between all the sequences on a MSA. You can indicate the output
element type with the last optional parameter (`Float64` by default). For a MSA with a lot
of sequences, you can use `Float32` or `Flot16` in order to avoid the `OutOfMemoryError()`.
"""
function percentidentity{T}(msa::AbstractMatrix{Residue}, out::Type{T}=Float64)
    aln = getresiduesequences(msa)
    nseq = length(aln)
    len = length(aln[1])
    scores = sequencepairsmatrix(msa, T, false, T(100.0))
    _percentidentity_kernel!(scores, aln, nseq, len)
    scores
end

# Mean Percent Identity of an MSA
# ===============================

"""
Returns the mean of the percent identity between the sequences of a MSA. If the MSA has 300
sequences or less, the mean is exact. If the MSA has more sequences and the `exact` keyword
is `false` (defualt), 44850 random pairs of sequences are used for the estimation. The
number of samples can be changed using the second argument. Use `exact=true` to perform all
the pairwise comparison (the calculation could be slow).
"""
function meanpercentidentity(msa, nsamples::Int=44850; exact::Bool=false)
    #                 lengthlist(300, false) == 44850
    nseq, ncol = size(msa)
    nvalues = lengthlist(nseq, Val{false})
    sum = 0.0
    if !exact && nvalues >= nsamples
        samples = Set{Tuple{Int,Int}}()
        sizehint!(samples, nsamples)
        while length(samples) < nsamples
            i = rand(1:(nseq-1))
            j = rand((i+1):nseq)
            if !in((i,j), samples)
                push!(samples, (i,j))
                sum += percentidentity(msa[i,:], msa[j,:])
            end
        end
        return( (sum/nsamples) )
    else # exact and/or few sequences
        for i in 1:(nseq-1)
            @inbounds @simd for j in (i+1):nseq
                sum += percentidentity(msa[i,:], msa[j,:])
            end
        end
        return( (sum/nvalues) )
    end
end

############################################################################################

# TO DO: Reduced alphabet

# Similarity percent
# ==================

"""
Calculates the similarity percent between two aligned sequences. The 100% is the length of
the aligned sequences minus the number of columns with gaps in both sequences.
Two residues are considered similar if they below to the same group.
The `groups` (third positional argument) can be indicated with a vector of length 20,
having the group labels of each `Residue` in the following order (mandatory):
A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V

**MIToS default groups** are:

```julia
[
:nonpolar, # A
:positive, # R
:polar,    # N
:negative, # D
:cysteine, # C
:polar,    # Q
:negative, # E
:glycine,  # G
:positive, # H
:nonpolar, # I
:nonpolar, # L
:positive, # K
:nonpolar, # M
:aromatic, # F
:proline,  # P
:polar,    # S
:polar,    # T
:aromatic, # W
:aromatic, # Y
:nonpolar  # V
]

```

Other residue groups:

**SMS (Sequence Manipulation Suite)** Ident and Sim: GAVLI, FYW, ST, KRH, DENQ, P, CM

*Stothard P (2000)
The Sequence Manipulation Suite: JavaScript programs for analyzing and formatting protein and DNA sequences.
Biotechniques 28:1102-1104.*

```julia

[1, 4, 5, 5, 7, 5, 5, 1, 4, 1, 1, 4, 7, 2, 6, 3, 3, 2, 2, 1]

```

**Bio3D 2.2** seqidentity: GA, MVLI, FYW, ST, KRH, DE, NQ, P, C

*Grant, B.J. et al. (2006) Bioinformatics 22, 2695--2696.*

```julia

[1, 5, 7, 6, 9, 7, 6, 1, 5, 2, 2, 5, 2, 3, 8, 4, 4, 3, 3, 2]

```
"""
function percentsimilarity(seq1::Vector{Residue}, seq2::Vector{Residue},
    alphabet::ResidueAlphabet = reduced"(AILMV)(RHK)(NQST)(DE)(FWY)CGP")

    len = length(seq1)
    if len != length(seq2)
        throw(ErrorException("""
        Sequences of different lengths, they aren't aligned or don't come from the same MSA.
        """))
    end

    count = 0
    colgap = 0
    colxaa = 0 # Columns with XAA or residues outside the alphabet aren't used
    @inbounds for i in 1:len
        res1 = seq1[i]
        res2 = seq2[i]
        if !in(res1, alphabet) || !in(res2, alphabet)
            colxaa += 1
            continue
        end
        isgap1 = res1 == GAP
        isgap2 = res2 == GAP
        if isgap1 && isgap2
            colgap += 1
        end
        if !isgap1 && !isgap2
            count += Int(alphabet[res1] == alphabet[res2])
        end
    end

    (100.0 * count)/(len-colgap-colxaa)
end

"""
Calculates the similarity percent between all the sequences on a MSA. You can indicate
the output element type with the `out` keyword argument (`Float64` by default). For an
MSA with a lot of sequences, you can use `out=Float32` or `out=Flot16` in order to
avoid the `OutOfMemoryError()`.
"""
function percentsimilarity(msa::AbstractMatrix{Residue}, A...; out::Type=Float64)
    M = getresiduesequences(msa)
    P = sequencepairsmatrix(msa, out, false, out(100.0))
    @inbounds @iterateupper P false list[k]=:($percentsimilarity)(:($M)[i],:($M)[j],:($A)...)
    P
end
