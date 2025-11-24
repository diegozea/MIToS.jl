"""
All is used instead of MIToS 1.0 "all" or "*", because it's possible to dispatch on it.
"""
struct All end

_get_function_name(str::AbstractString)::String = split(str, '.')[end]

"""
`get_n_words{T <: Union{ASCIIString, UTF8String}}(line::T, n::Int)`
It returns a `Vector{T}` with the first `n` (possibles) words/fields (delimited
by space or tab). If there is more than `n` words, the last word
returned contains the finals words and the delimiters. The length of the
returned vector is `n` or less (if the number of words is less than `n`).
This is used for parsing the Stockholm format.

```jldoctest
julia> using MIToS.Utils

julia> get_n_words("#=GR O31698/18-71 SS    CCCHHHHHHHHHHHHHHHEEEEEEEEEEEEEEEEHHH", 3)
3-element Vector{String}:
 "#=GR"
 "O31698/18-71"
 "SS    CCCHHHHHHHHHHHHHHHEEEEEEEEEEEEEEEEHHH"
```
"""
function get_n_words(line::String, n::Int)
    if isempty(line)
        return String[]
    end
    idx = firstindex(line)
    last_idx = lastindex(line)

    # Skip leading delimiters so we do not emit empty fields for padded lines.
    while idx <= last_idx
        c = @inbounds line[idx]
        if c == ' ' || c == '\t'
            idx = nextind(line, idx)
        else
            break
        end
    end

    if idx > last_idx
        return String[]
    end

    words = Array{String}(undef, n)
    start_idx = idx
    prev_idx = idx  # remember the last character boundary to avoid repeated prevind calls
    word_count = 0
    max_split = n - 1

    # Collect at most the first n - 1 words; the nth entry keeps the remainder.
    while word_count < max_split && idx <= last_idx
        c = @inbounds line[idx]
        if c == ' ' || c == '\t'
            word_count += 1
            @inbounds words[word_count] = line[start_idx:prev_idx]

            # Skip the complete delimiter run.
            idx = nextind(line, idx)
            while idx <= last_idx
                c = @inbounds line[idx]
                if c == ' ' || c == '\t'
                    idx = nextind(line, idx)
                else
                    break
                end
            end
            start_idx = idx
            prev_idx = idx
            continue
        end
        prev_idx = idx
        idx = nextind(line, idx)
    end

    if start_idx <= last_idx && word_count < n
        word_count += 1
        @inbounds words[word_count] = line[start_idx:last_idx]
    end

    word_count < n && resize!(words, word_count)
    return words
end

"""
`hascoordinates(id)`
It returns `true` if `id`/sequence name has the format: **UniProt/start-end**
(i.e. O83071/192-246)
"""
function hascoordinates(id)
    occursin(r"^\w+/\d+-\d+$", id)
end

"""
Selects the first element of the vector. This is useful for unpacking one element vectors.
Throws a warning if there are more elements. `element_name` is *element* by default,
but the name can be changed using the second argument.
"""
function select_element(
    vector::Array{T,1},
    element_name::AbstractString = "element",
) where {T}
    len = length(vector)
    if len == 0
        throw(ErrorException("There is not $element_name"))
    elseif len != 1
        @warn("There are more than one ($len) $element_name using the first.")
    end
    @inbounds return (vector[1])
end

"""
Returns a vector with the `part` ("upper" or "lower") of the square matrix `mat`.
The `diagonal` is not included by default.
"""
function matrix2list(
    mat::AbstractMatrix{T};
    part = "upper",
    diagonal::Bool = false,
) where {T}
    nrow, ncol = size(mat)
    if nrow != ncol
        throw(ErrorException("Should be a square matrix"))
    end
    if diagonal
        d = 0
        N = div((ncol * ncol) + ncol, 2)
    else
        d = 1
        N = div((ncol * ncol) - ncol, 2)
    end
    list = Array{T}(undef, N)
    k = 1
    if part == "upper"
        for i = 1:(ncol-d)
            for j = (i+d):ncol
                list[k] = mat[i, j]
                k += 1
            end
        end
    elseif part == "lower"
        for j = 1:(ncol-d)
            for i = (j+d):ncol
                list[k] = mat[i, j]
                k += 1
            end
        end
    else
        throw(ErrorException("part should be \"upper\" or \"lower\""))
    end
    list
end

"""
Returns a square symmetric matrix from the vector `vec`. `side` is the number of
rows/columns. The `diagonal` is not included by default, set to `true` if there are
diagonal elements in the list.
"""
function list2matrix(vec::AbstractVector{T}, side::Int; diagonal::Bool = false) where {T}
    d = diagonal ? 0 : 1
    mat = zeros(T, side, side)
    k = 1
    for i = 1:(side-d)
        for j = (i+d):side
            value = vec[k]
            mat[i, j] = value
            mat[j, i] = value
            k += 1
        end
    end
    mat
end

"""
It checks if a PDB code has the correct format.
"""
check_pdbcode(pdbcode::AbstractString) = occursin(r"^\w{4}$", pdbcode)

"""
Getter for the `array` field of `NamedArray`s
"""
getarray(x::NamedArray) = x.array
