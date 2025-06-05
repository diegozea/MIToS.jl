const _BLOSUM62_DATA = Int[
    #    A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   -   X
    4  -1  -2  -2   0  -1  -1   0  -2  -1  -1  -1  -1  -2  -1   1   0  -3  -2   0  -4  -1
   -1   5   0  -2  -3   1   0  -2   0  -3  -2   2  -1  -3  -2  -1  -1  -3  -2  -3  -4  -1
   -2   0   6   1  -3   0   0   0   1  -3  -3   0  -2  -3  -2   1   0  -4  -2  -3  -4  -1
   -2  -2   1   6  -3   0   2  -1  -1  -3  -4  -1  -3  -3  -1   0  -1  -4  -3  -3  -4  -1
    0  -3  -3  -3   9  -3  -4  -3  -3  -1  -1  -3  -1  -2  -3  -1  -1  -2  -2  -1  -4  -1
   -1   1   0   0  -3   5   2  -2   0  -3  -2   1   0  -3  -1   0  -1  -2  -1  -2  -4  -1
   -1   0   0   2  -4   2   5  -2   0  -3  -3   1  -2  -3  -1   0  -1  -3  -2  -2  -4  -1
    0  -2   0  -1  -3  -2  -2   6  -2  -4  -4  -2  -3  -3  -2   0  -2  -2  -3  -3  -4  -1
   -2   0   1  -1  -3   0   0  -2   8  -3  -3  -1  -2  -1  -2  -1  -2  -2   2  -3  -4  -1
   -1  -3  -3  -3  -1  -3  -3  -4  -3   4   2  -3   1   0  -3  -2  -1  -3  -1   3  -4  -1
   -1  -2  -3  -4  -1  -2  -3  -4  -3   2   4  -2   2   0  -3  -2  -1  -2  -1   1  -4  -1
   -1   2   0  -1  -3   1   1  -2  -1  -3  -2   5  -1  -3  -1   0  -1  -3  -2  -2  -4  -1
   -1  -1  -2  -3  -1   0  -2  -3  -2   1   2  -1   5   0  -2  -1  -1  -1  -1   1  -4  -1
   -2  -3  -3  -3  -2  -3  -3  -3  -1   0   0  -3   0   6  -4  -2  -2   1   3  -1  -4  -1
   -1  -2  -2  -1  -3  -1  -1  -2  -2  -3  -3  -1  -2  -4   7  -1  -1  -4  -3  -2  -4  -1
    1  -1   1   0  -1   0   0   0  -1  -2  -2   0  -1  -2  -1   4   1  -3  -2  -2  -4  -1
    0  -1   0  -1  -1  -1  -1  -2  -2  -1  -1  -1  -1  -2  -1   1   5  -2  -2   0  -4  -1
   -3  -3  -4  -4  -2  -2  -3  -2  -2  -3  -2  -3  -1   1  -4  -3  -2  11   2  -3  -4  -1
   -2  -2  -2  -3  -2  -1  -2  -3   2  -1  -1  -2  -1   3  -3  -2  -2   2   7  -1  -4  -1
    0  -3  -3  -3  -1  -2  -2  -3  -3   3   1  -2   1  -1  -2  -2   0  -3  -1   4  -4  -1
   -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4   1  -4
   -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -4  -1
] # 22x22

const BLOSUM62_SCORE = _BLOSUM62_DATA

function _pairwise_score(
    seq1::AbstractVector{Residue},
    seq2::AbstractVector{Residue};
    gap_open::Int,
    gap_extend::Int,
    matrix::AbstractMatrix{Int} = BLOSUM62_SCORE,
)
    len = length(seq1)
    score = 0
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
            score += matrix[Int(r1), Int(r2)]
        end
    end
    score
end

"""
    sum_of_pairs_score(msa; gap_open=-10, gap_extend=-1, matrix=BLOSUM62_SCORE)

Calculate the sum-of-pairs score of an alignment using a substitution matrix.
By default, the BLOSUM62 matrix is used. The penalties for gap opening and
extension can be adjusted with the `gap_open` and `gap_extend` keyword
arguments. Pass a different `matrix` to use another substitution scheme.
"""
function sum_of_pairs_score(
    msa::AbstractMatrix{Residue};
    gap_open::Int = -10,
    gap_extend::Int = -1,
    matrix::AbstractMatrix{Int} = BLOSUM62_SCORE,
)
    nseq, _ = size(msa)
    total = 0
    @inbounds for i = 1:(nseq-1)
        seqi = view(msa, i, :)
        for j = (i+1):nseq
            total += _pairwise_score(
                seqi,
                view(msa, j, :);
                gap_open = gap_open,
                gap_extend = gap_extend,
                matrix = matrix,
            )
        end
    end
    total
end

sum_of_pairs_score(msa::AbstractAlignedObject; kargs...) =
    sum_of_pairs_score(namedmatrix(msa); kargs...)
