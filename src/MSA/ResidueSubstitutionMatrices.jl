module ResidueSubstitutionMatrices

using ..MSA: Residue, ResidueAlphabet, GappedAlphabet
using Printf

export ResidueSubstitutionMatrix

"""
    ResidueSubstitutionMatrix{T<:Real}(scores::Matrix{T}, alphabet::ResidueAlphabet = GappedAlphabet())

Symmetric substitution scores for pairs of [`Residue`](@ref) values.
The matrix must be square and its size must match the length of the
provided alphabet. By default a [`GappedAlphabet`](@ref) is used (the
`'*'` character of classical substitution matrices is represented by `'-'`).
"""
struct ResidueSubstitutionMatrix{T<:Real,A<:ResidueAlphabet} <: AbstractMatrix{T}
    scores::Matrix{T}
    alphabet::A

    function (::Type{ResidueSubstitutionMatrix{T,A}})(
        scores::Matrix{T},
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

ResidueSubstitutionMatrix(scores::Matrix{T}) where {T<:Real} =
    ResidueSubstitutionMatrix(scores, GappedAlphabet())

ResidueSubstitutionMatrix(
    scores::Matrix{T},
    alphabet::A,
) where {T<:Real,A<:ResidueAlphabet} = ResidueSubstitutionMatrix{T,A}(scores, alphabet)

"""
    matrix[a::Residue, b::Residue]

Return the substitution score for the pair of residues `a` and `b`.
"""
@inline function Base.getindex(
@inline function Base.getindex(
    matrix::ResidueSubstitutionMatrix{T},
    r::Residue,
    c::Colon,
) where {T}
    i = matrix.alphabet[r]
    @inbounds matrix.scores[i, c]
end

@inline function Base.getindex(
    matrix::ResidueSubstitutionMatrix{T},
    r::Colon,
    c::Residue,
) where {T}
    j = matrix.alphabet[c]
    @inbounds matrix.scores[r, j]
end

    labels = names(matrix.alphabet)
    label_width = maximum(length.(labels))
    value_strings = [@sprintf("%g", v) for v in matrix.scores]
    value_width = maximum(length.(value_strings)) + 1
    n = size(matrix.scores, 1)
    println(io, typeof(matrix), " :")
    for i = 1:n
        print(io, rpad(labels[i], label_width), " ")
        offset = (i - 1) * n
        for j = 1:n
            idx = offset + j
            print(io, lpad(value_strings[idx], value_width))
        end
        i == n || println(io)
    end
    i = matrix.alphabet[a]
    j = matrix.alphabet[b]
    @inbounds matrix.scores[i, j]
end

function Base.show(io::IO, ::MIME"text/plain", matrix::ResidueSubstitutionMatrix)
    println(io, typeof(matrix), " : ")
    println(io, "Alphabet: ", matrix.alphabet)
    show(io, matrix.scores)
end

end # module
