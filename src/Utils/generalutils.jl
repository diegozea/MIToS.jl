"""
`deleteitems!(vector, items)` deletes a group of `items` (test it using `in`) from `vector`.
"""
function deleteitems!(vector::Vector, items)
  i = 1
  j = 0
  @inbounds while j < length(vector)
    j += 1
    value = vector[j]
    if !(value in items)
      if i != j
        vector[i] = value
      end
      i += 1
    end
  end
  resize!(vector, i-1)
end

_get_data(str::ASCIIString) = str.data
_get_data(str::UTF8String) = collect(str)

"""
`get_n_words{T <: Union{ASCIIString, UTF8String}}(line::T, n::Int)` returns a `Vector{T}` with the first `n` (possibles) words/fields (delimited by space, tab or newline).
If there is more than `n` words, the last word returned contains the finals words and the delimiters.
The length of the returned vector is `n` or less (if the number of words is less than `n`).
The last newline character is always removed.
This is used for parsing the Stockholm format.
```
julia> get_n_words("#=GR O31698/18-71 SS    CCCHHHHHHHHHHHHHHHEEEEEEEEEEEEEEEEHHH\n", 3)
3-element Array{ASCIIString,1}:
 "#=GR"
 "O31698/18-71"
 "SS    CCCHHHHHHHHHHHHHHHEEEEEEEEEEEEEEEEHHH"
```
"""
function get_n_words{T <: ByteString}(line::T, n::Int)
  linedata = _get_data(line)
  CharType = eltype(linedata)
  N = length(linedata)
  coords = Array(T, n)
  spaces = CharType[' ', '\n', '\t']
  j = 0
  start = 0
  ended = false
  if N > 1
    block = !(linedata[1] in spaces)
    if block
      start = 1
    end
    for i in 2:N-1
      if linedata[i] in spaces
        if block
          j += 1
          if j == n
            last = linedata[N] == CharType('\n') ? N-1 : N
            @inbounds coords[j] = UTF8String(getindex(linedata,start:last))
            ended = true
            break
          else
            @inbounds coords[j] = UTF8String(getindex(linedata,start:i-1))
          end
        end
        block = false
      else
        if ! block
          block = true
          start = i
        end
      end
    end
    if !ended
      j += 1
      last = linedata[N] == CharType('\n') ? N-1 : N
      @inbounds coords[j] = UTF8String(getindex(linedata,start:last))
    end
  elseif N==1
    if linedata[1] in spaces
      j=0
    else
      j=1
      coords[1] = line
    end
  else
    j=0
  end
  if j != n
    resize!(coords, j)
  end
  coords
end

"""
`hascoordinates(id)` returns `true` if `id` has the format: **uniprot/start-end** i.e. O83071/192-246
"""
function hascoordinates(id)
  ismatch(r"^\w+/\d+-\d+$", id)
end

"""
Selects the first element of the vector.
This is useful for unpacking one element vectors.
Throws a warning if there are more elements.
The `element_name` for the warning is *element* by default, but the name can be changed.
"""
function select_element{T}(vector::Array{T,1}, element_name::ASCIIString="element")
  len = length(vector)
  if len == 0
    throw(ErrorException(string("There is not ", element_name)))
  elseif len != 1
    warn(string("There are more than one (", len, ") ", element_name, " using the first."))
  end
  @inbounds return(vector[1])
end

"""
Returns a vector with the `part` ("upper" or "lower") of the square matrix `mat`.
The `diagonal` is not included by default.
"""
function matrix2list{T}(mat::AbstractMatrix{T}; part="upper", diagonal::Bool=false)
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
  list = Array(T, N)
  k = 1
  if part=="upper"
    for i in 1:(ncol-d)
      for j in (i+d):ncol
        list[k] = mat[i, j]
        k += 1
      end
    end
  elseif part=="lower"
    for j in 1:(ncol-d)
      for i in (j+d):ncol
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
Returns a square symmetric matrix from the vector `vec`.
`side` is the number of rows/columns.
The `diagonal` is not included by default, set to `true` if the diagonal elements are in the list.
"""
function list2matrix{T}(vec::AbstractVector{T}, side::Int; diagonal::Bool=false)
  d = diagonal ? 0 : 1
  mat = zeros(T, side, side)
  k = 1
  for i in 1:(side-d)
    for j in (i+d):side
      value = vec[k]
      mat[i, j] = value
      mat[j, i] = value
      k += 1
    end
  end
  mat
end
