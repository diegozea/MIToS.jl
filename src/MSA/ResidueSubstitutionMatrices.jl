module ResidueSubstitutionMatrices

using ..MSA: Residue, ResidueAlphabet, GappedAlphabet
using Printf

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
    labels = names(matrix.alphabet)
    strlabels = map(labels) do lbl
        length(lbl) == 1 ? lbl : "(" * lbl * ")"
    end
    width = maximum(length.(strlabels))
    for (i, lbl) in enumerate(strlabels)
        print(io, rpad(lbl, width), " ")
        row = matrix.scores[i, :]
        for v in row
            Printf.@printf(io, "%8g", v)
        end
        i < length(strlabels) && print(io, '\n')
    end
end

Base.show(io::IO, matrix::ResidueSubstitutionMatrix) = show(io, MIME"text/plain"(), matrix)

end # module