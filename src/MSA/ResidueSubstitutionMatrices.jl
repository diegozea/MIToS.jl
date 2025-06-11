module ResidueSubstitutionMatrices

using ..MSA: Residue, ResidueAlphabet, GappedAlphabet, getnamedict
using NamedArrays

export ResidueSubstitutionMatrix

"""
    ResidueSubstitutionMatrix{T<:Real}(scores::AbstractMatrix{T}, alphabet::ResidueAlphabet = GappedAlphabet())

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

ResidueSubstitutionMatrix(scores::AbstractMatrix{T}) where {T<:Real} =
    ResidueSubstitutionMatrix(scores, GappedAlphabet())

ResidueSubstitutionMatrix(
    scores::AbstractMatrix{T},
    alphabet::A,
) where {T<:Real,A<:ResidueAlphabet} = ResidueSubstitutionMatrix{T,A}(scores, alphabet)

"""
    matrix[i::Int, j::Int]

Return the substitution score for the pair of residue indexes `i` and `j`.
"""
@inline Base.getindex(matrix::ResidueSubstitutionMatrix{T}, i::Int, j::Int) where {T} =
    matrix.scores[i, j]

"""
    matrix[a::Residue, b::Residue]

Return the substitution score for the pair of residues `a` and `b`.
"""
@inline function Base.getindex(
    matrix::ResidueSubstitutionMatrix{T},
    a::Residue,
    b::Residue,
) where {T}
    i = matrix.alphabet[a]
    j = matrix.alphabet[b]
    @inbounds matrix.scores[i, j]
end

function Base.show(io::IO, ::MIME"text/plain", matrix::ResidueSubstitutionMatrix)
    labels = getnamedict(matrix.alphabet)
    named = NamedArray(matrix.scores, (labels, labels), ("Residue", "Residue"))
    print(io, typeof(matrix), " : ")
    show(io, MIME"text/plain"(), named)
end

end # module
