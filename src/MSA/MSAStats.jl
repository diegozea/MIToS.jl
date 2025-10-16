
# Counting Gaps and Coverage
# --------------------------

"""
It calculates the fraction of gaps on the `Array` (alignment, sequence, column, etc.). This
function can take an extra `dimension` argument for calculation of the gap fraction over the
given dimension.
"""
function gapfraction(x::AbstractArray{Residue})
    counter = 0
    len = 0
    for res in x
        counter += ifelse(res == GAP, 1, 0)
        len += 1
    end
    float(counter) / float(len)
end

function gapfraction(x::AbstractArray{Residue}, dimension::Int)
    mapslices(gapfraction, x, dims = dimension)
end

"""
It calculates the fraction of residues (no gaps) on the `Array`
(alignment, sequence, column, etc.). This function can take an extra `dimension` argument
for calculation of the residue fraction over the given dimension
"""
residuefraction(x::AbstractArray{Residue}) = 1.0 - gapfraction(x)

function residuefraction(x::AbstractArray{Residue}, dimension::Int)
    mapslices(residuefraction, x, dims = dimension)
end

macro keep_names_dimension(functions)
    function_names = functions.args
    n = length(function_names)
    definitions = Array{Any}(undef, n)

    for i = 1:n
        f = esc(function_names[i])
        definitions[i] = quote

            function ($f)(msa::NamedResidueMatrix{T}, dimension::Int) where {T}
                @assert dimension == 1 || dimension == 2 "Dimension must be 1 or 2."
                result = ($f)(getarray(msa), dimension)
                if dimension == 1
                    name_list = names(msa, 2)
                    N = length(name_list)
                    NamedArray(
                        result,
                        (
                            OrderedDict{String,Int}(
                                Utils._get_function_name(string($f)) => 1,
                            ),
                            OrderedDict{String,Int}(name_list[i] => i for i = 1:N),
                        ),
                        ("Function", "Col"),
                    )
                elseif dimension == 2
                    name_list = names(msa, 1)
                    N = length(name_list)
                    NamedArray(
                        result,
                        (
                            OrderedDict{String,Int}(name_list[i] => i for i = 1:N),
                            OrderedDict{String,Int}(
                                Utils._get_function_name(string($f)) => 1,
                            ),
                        ),
                        ("Seq", "Function"),
                    )
                end
            end

            ($f)(a::AbstractResidueMatrix, dimension::Int) =
                ($f)(namedmatrix(a), dimension)
        end
    end

    return Expr(:block, definitions...)
end

@keep_names_dimension([gapfraction, residuefraction])

"""
Coverage of the sequences with respect of the number of positions on the MSA
"""
function coverage(msa::AbstractMatrix{Residue})
    result = residuefraction(msa, 2)
    if isa(result, NamedArray) && ndims(result) == 2
        setnames!(result, ["coverage"], 2)
    end
    result
end

coverage(msa::AbstractAlignedObject) = coverage(namedmatrix(msa))

"""
Fraction of gaps per column/position on the MSA
"""
columngapfraction(msa::AbstractMatrix{Residue}) = gapfraction(msa, 1)
columngapfraction(msa::AbstractAlignedObject) = columngapfraction(namedmatrix(msa))

"""
    n_effective(msa[, threshold=80.0])

Calculate the effective number of sequences (`Neff`/`Meff`) of an alignment using
the weighting scheme described by *Morcos et al.*.
For each sequence the weight is the inverse of the number of sequences
(including itself) whose percent identity is greater than or equal to
`threshold`. The sum of the weights is returned.

# References

  - [Morcos, Faruck, et al. "Direct-coupling analysis of residue coevolution
    captures native contacts across many protein families." Proceedings of the
    National Academy of Sciences 108.49 (2011): E1293–E1301.](@cite 10.1073/pnas.1111471108)
"""
function n_effective(msa::AbstractMatrix{Residue}, threshold::Real = 80.0)
    sequences = getresiduesequences(msa)
    nseq = length(sequences)
    P = sequencepairsmatrix(msa, Bool, Val{false}, true)
    @inbounds @iterateupper getarray(P) false begin
        list[k] = percentidentity(sequences[i], sequences[j], threshold)
    end
    counts = fill(1, nseq)
    @inbounds @iterateupper getarray(P) false begin
        if list[k]
            counts[i] += 1
            counts[j] += 1
        end
    end
    total = 0.0
    @inbounds for c in counts
        total += 1.0 / c
    end
    total
end

n_effective(msa::AbstractAlignedObject, threshold::Real = 80.0) =
    n_effective(namedmatrix(msa), threshold)
