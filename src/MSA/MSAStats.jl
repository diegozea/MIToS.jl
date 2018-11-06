
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
    mapslices(gapfraction, x, dims=dimension)
end

"""
It calculates the fraction of residues (no gaps) on the `Array`
(alignment, sequence, column, etc.). This function can take an extra `dimension` argument
for calculation of the residue fraction over the given dimension
"""
residuefraction(x::AbstractArray{Residue}) = 1.0 - gapfraction(x)

function residuefraction(x::AbstractArray{Residue}, dimension::Int)
    mapslices(residuefraction, x, dims=dimension)
end

macro keep_names_dimension(functions)
    function_names = functions.args
    n = length(function_names)
    definitions = Array{Any}(undef, n)

    for i in 1:n
        f = esc(function_names[i])
        definitions[i] = quote

            function ($f)(msa::NamedResidueMatrix{T}, dimension::Int) where T
                result = ($f)(getarray(msa), dimension)
                if dimension == 1
                    name_list = names(msa,2)
                    N = length(name_list)
                    NamedArray(result,
                        (OrderedDict{String,Int}(Utils._get_function_name(string($f))=>1),
                         OrderedDict{String,Int}(name_list[i]=>i for i in 1:N)),
                        ("Function","Col"))
                elseif dimension == 2
                    name_list = names(msa,1)
                    N = length(name_list)
                    NamedArray(result,
                        (OrderedDict{String,Int}(name_list[i]=>i for i in 1:N),
                         OrderedDict{String,Int}(Utils._get_function_name(string($f))=>1)),
                        ("Seq","Function"))
                else
                    throw(ArgumentError("Dimension must be 1 or 2."))
                end
            end

            ($f)(a::AbstractAlignedObject, dimension::Int) = ($f)(namedmatrix(a), dimension)
        end
    end

    return Expr(:block, definitions...)
end

@keep_names_dimension([gapfraction, residuefraction])

"Coverage of the sequences with respect of the number of positions on the MSA"
function coverage(msa::AbstractMatrix{Residue})
    result = residuefraction(msa, 2)
    if isa(result, NamedArray) && ndims(result) == 2
        setnames!(result,["coverage"],2)
    end
    result
end

coverage(msa::AbstractAlignedObject) = coverage(namedmatrix(msa))

"Fraction of gaps per column/position on the MSA"
columngapfraction(msa::AbstractMatrix{Residue}) = gapfraction(msa, 1)
columngapfraction(msa::AbstractAlignedObject) = columngapfraction(namedmatrix(msa))
