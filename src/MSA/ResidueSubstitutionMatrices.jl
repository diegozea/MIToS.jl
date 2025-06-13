module ResidueSubstitutionMatrices

using ..MSA: Residue, ResidueAlphabet, GappedAlphabet, GappedXAlphabet
import ..MSA: getnamedict
using NamedArrays

export ResidueSubstitutionMatrix


"""
    ResidueSubstitutionMatrix{T<:Real, A<:ResidueAlphabet}(
        scores::AbstractMatrix{T},
        alphabet::A = GappedAlphabet(),
    )

Symmetric substitution scores for pairs of [`Residue`](@ref) values.
The matrix must be square and its size must match the length of the
provided alphabet. By default a [`GappedAlphabet`](@ref) is used (the
`'*'` character of classical substitution matrices is represented by `'-'`).
"""
struct ResidueSubstitutionMatrix{T<:Real,A<:ResidueAlphabet} <: AbstractMatrix{T}
    scores::Matrix{T}
    alphabet::A

    function (::Type{ResidueSubstitutionMatrix{T,A}})(
        scores::AbstractMatrix{T},
        alphabet::A,
    ) where {T<:Real,A<:ResidueAlphabet}
        n, m = size(scores)
        if n != m
            throw(ArgumentError("scores matrix must be square"))
        end
        if n != length(alphabet)
            throw(ArgumentError("matrix size must match alphabet length"))
        end
        new{T,A}(Matrix{T}(scores), alphabet)
    end
end

Base.size(matrix::ResidueSubstitutionMatrix) = size(matrix.scores)

@inline Base.getindex(matrix::ResidueSubstitutionMatrix, I::Vararg{Int,N}) where {N} =
    matrix.scores[I...]

ResidueSubstitutionMatrix(
    scores::AbstractMatrix{T},
    alphabet::A = GappedAlphabet(),
) where {T<:Real,A<:ResidueAlphabet} = ResidueSubstitutionMatrix{T,A}(scores, alphabet)


"""
    matrix[a::Residue, b::Residue]

Return the substitution score for the pair of residues `a` and `b`.
"""
@inline function Base.getindex(
    matrix::ResidueSubstitutionMatrix{T,A},
    a::Residue,
    b::Residue,
) where {T,A}
    i = matrix.alphabet[a]
    j = matrix.alphabet[b]
    @inbounds matrix.scores[i, j]
end

@inline function Base.getindex(
    matrix::ResidueSubstitutionMatrix{T,A},
    r::Residue,
    c::Colon,
) where {T,A}
    i = matrix.alphabet[r]
    @inbounds matrix.scores[i, c]
end

@inline function Base.getindex(
    matrix::ResidueSubstitutionMatrix{T,A},
    r::Colon,
    c::Residue,
) where {T,A}
    j = matrix.alphabet[c]
    @inbounds matrix.scores[r, j]
end

function Base.show(io::IO, ::MIME"text/plain", matrix::ResidueSubstitutionMatrix)
    println(io, "$(typeof(matrix)) :")
    rawlabels = names(matrix.alphabet)
    labels = [length(lbl) == 1 ? lbl : "(" * lbl * ")" for lbl in rawlabels]
    named = NamedArray(matrix.scores, (labels, labels))
    ctx = IOContext(io, :limit => true)
    show(ctx, MIME"text/plain"(), named, size(named, 1))
end

Base.show(io::IO, matrix::ResidueSubstitutionMatrix) = show(io, MIME"text/plain"(), matrix)

# This BLOSUM62 matrix was obtained from the NCBI source code :
# https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/BLOSUM62
"""
    BLOSUM62

Standard BLOSUM62 substitution matrix represented as a
[`ResidueSubstitutionMatrix`](@ref) using [`GappedXAlphabet`](@ref).
"""
const BLOSUM62 = ResidueSubstitutionMatrix(
    Int[
        #    A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T
        #    W   Y   V   -   X
        4 -1 -2 -2 0 -1 -1 0 -2 -1 -1 -1 -1 -2 -1 1 0 -3 -2 0 -4 -1
        -1 5 0 -2 -3 1 0 -2 0 -3 -2 2 -1 -3 -2 -1 -1 -3 -2 -3 -4 -1
        -2 0 6 1 -3 0 0 0 1 -3 -3 0 -2 -3 -2 1 0 -4 -2 -3 -4 -1
        -2 -2 1 6 -3 0 2 -1 -1 -3 -4 -1 -3 -3 -1 0 -1 -4 -3 -3 -4 -1
        0 -3 -3 -3 9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -4 -1
        -1 1 0 0 -3 5 2 -2 0 -3 -2 1 0 -3 -1 0 -1 -2 -1 -2 -4 -1
        -1 0 0 2 -4 2 5 -2 0 -3 -3 1 -2 -3 -1 0 -1 -3 -2 -2 -4 -1
        0 -2 0 -1 -3 -2 -2 6 -2 -4 -4 -2 -3 -3 -2 0 -2 -2 -3 -3 -4 -1
        -2 0 1 -1 -3 0 0 -2 8 -3 -3 -1 -2 -1 -2 -1 -2 -2 2 -3 -4 -1
        -1 -3 -3 -3 -1 -3 -3 -4 -3 4 2 -3 1 0 -3 -2 -1 -3 -1 3 -4 -1
        -1 -2 -3 -4 -1 -2 -3 -4 -3 2 4 -2 2 0 -3 -2 -1 -2 -1 1 -4 -1
        -1 2 0 -1 -3 1 1 -2 -1 -3 -2 5 -1 -3 -1 0 -1 -3 -2 -2 -4 -1
        -1 -1 -2 -3 -1 0 -2 -3 -2 1 2 -1 5 0 -2 -1 -1 -1 -1 1 -4 -1
        -2 -3 -3 -3 -2 -3 -3 -3 -1 0 0 -3 0 6 -4 -2 -2 1 3 -1 -4 -1
        -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4 7 -1 -1 -4 -3 -2 -4 -1
        1 -1 1 0 -1 0 0 0 -1 -2 -2 0 -1 -2 -1 4 1 -3 -2 -2 -4 -1
        0 -1 0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1 1 5 -2 -2 0 -4 -1
        -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1 1 -4 -3 -2 11 2 -3 -4 -1
        -2 -2 -2 -3 -2 -1 -2 -3 2 -1 -1 -2 -1 3 -3 -2 -2 2 7 -1 -4 -1
        0 -3 -3 -3 -1 -2 -2 -3 -3 3 1 -2 1 -1 -2 -2 0 -3 -1 4 -4 -1
        -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 1 -4
        -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -4 -1
    ],
    GappedXAlphabet(),
)

end # module
