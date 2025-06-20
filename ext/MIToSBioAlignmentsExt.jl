module MIToSBioAlignmentsExt

using MIToS.MSA
using MIToS.MSA.ResidueSubstitutionMatrices: ResidueSubstitutionMatrix

import BioAlignments

function Base.convert(
    ::Type{ResidueSubstitutionMatrix{T,A}},
    submat::BioAlignments.SubstitutionMatrix{BioAlignments.BioSymbols.AminoAcid,S},
) where {T<:Real,A<:ResidueAlphabet,S<:Real}
    alphabet = A()
    n = length(alphabet)
    scores = fill(zero(T), n, n)
    aa_alphabet = BioAlignments.BioSymbols.alphabet(submat)
    for aa in aa_alphabet
        c_a = Char(aa)
        i = alphabet[Residue(c_a)]
        for bb in aa_alphabet
            c_b = Char(bb)
            j = alphabet[Residue(c_b)]
            scores[i, j] = convert(T, submat[c_a, c_b])
        end
    end
    gap = alphabet[GAP]
    scores[gap, :] .= -4
    scores[:, gap] .= -4
    scores[gap, gap] = 1

    x = alphabet[XAA]
    scores[x, :] .= -1
    scores[:, x] .= -1
    scores[x, x] = -1
    scores[gap, x] = -4
    scores[x, gap] = -4
    return ResidueSubstitutionMatrix{T,A}(scores, alphabet)
end

Base.convert(
    ::Type{ResidueSubstitutionMatrix},
    submat::BioAlignments.SubstitutionMatrix{BioAlignments.BioSymbols.AminoAcid,S},
) where {S<:Real} = convert(ResidueSubstitutionMatrix{S,GappedXAlphabet}, submat)

end
