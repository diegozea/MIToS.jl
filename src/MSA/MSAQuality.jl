
function _pairwise_score(
    seq1::AbstractVector{Residue},
    seq2::AbstractVector{Residue},
    gap_open::Real,
    gap_extend::Real,
    matrix::ResidueSubstitutionMatrices.ResidueSubstitutionMatrix,
)
    score = zero(typeof(gap_open))
    len = length(seq1)
    ingap1 = false
    ingap2 = false
    @inbounds for k = 1:len
        r1 = seq1[k]
        r2 = seq2[k]
        if r1 == GAP && r2 == GAP
            continue
        elseif r1 == GAP
            if ingap1
                score += gap_extend
            else
                score += gap_open
                ingap1 = true
            end
            ingap2 = false
        elseif r2 == GAP
            if ingap2
                score += gap_extend
            else
                score += gap_open
                ingap2 = true
            end
            ingap1 = false
        else
            ingap1 = false
            ingap2 = false
            score += convert(typeof(gap_open), matrix[r1, r2])
        end
    end
    score
end

"""
    sum_of_pairs_score(msa; gap_open=-10, gap_extend=-1,
                       matrix=ResidueSubstitutionMatrices.BLOSUM62)

Calculate the sum-of-pairs score of an alignment using a substitution matrix.
By default, the BLOSUM62 matrix is used. The penalties for gap opening and
extension can be adjusted with the `gap_open` and `gap_extend` keyword
arguments. Pass a different `matrix` to use another substitution scheme.
"""
function sum_of_pairs_score(
    msa::AbstractMatrix{Residue};
    gap_open::Real = -10,
    gap_extend::Real = -1,
    matrix::ResidueSubstitutionMatrices.ResidueSubstitutionMatrix{T,A} = ResidueSubstitutionMatrices.BLOSUM62,
) where {T,A}
    S = promote_type(T, typeof(gap_open), typeof(gap_extend))
    go = convert(S, gap_open)
    ge = convert(S, gap_extend)
    nseq, _ = size(msa)
    total = zero(S)
    @inbounds for i = 1:(nseq-1)
        seqi = view(msa, i, :)
        for j = (i+1):nseq
            total += _pairwise_score(seqi, view(msa, j, :), go, ge, matrix)
        end
    end
    total
end

sum_of_pairs_score(msa::AbstractAlignedObject; kargs...) =
    sum_of_pairs_score(namedmatrix(msa); kargs...)
