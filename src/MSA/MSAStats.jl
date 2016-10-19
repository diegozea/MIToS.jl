
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
    mapslices(gapfraction, x, dimension)
end

"""
It calculates the fraction of residues (no gaps) on the `Array`
(alignment, sequence, column, etc.). This function can take an extra `dimension` argument
for calculation of the residue fraction over the given dimension
"""
residuefraction(x::AbstractArray{Residue}) = 1.0 - gapfraction(x)

function residuefraction(x::AbstractArray{Residue}, dimension::Int)
    mapslices(residuefraction, x, dimension)
end

# The next functions keeps column or sequence names:
for f in (:gapfraction, :residuefraction)
    @eval begin
        function ($f)(msa::NamedArray{Residue}, dimension::Int)
            result = ($f)(array(msa), dimension)
            if dimension == 1
                NamedArray( result,
                            (OrderedDict([string($f)]), OrderedDict(names(msa,2))),
                            ("Function","Col"))
            elseif dimension == 2
                NamedArray( result,
                            (OrderedDict([string($f)]), OrderedDict(names(msa,1))),
                            ("Seq","Function"))
            else
                throw(ArgumentError("Dimension must be 1 or 2."))
            end
        end
        @eval ($f)(a::AbstractAlignedObject, dimension::Int) = ($f)(namedmatrix(a), d)
    end
end

"Coverage of the sequences with respect of the number of positions on the MSA"
coverage(msa::AbstractMatrix{Residue}) = residuefraction(msa, 2)
coverage(msa::AbstractAlignedObject) = coverage(namedmatrix(msa))

"Fraction of gaps per column/position on the MSA"
columngapfraction(msa::AbstractMatrix{Residue}) = gapfraction(msa, 1)
columngapfraction(msa::AbstractAlignedObject) = columngapfraction(namedmatrix(msa))
