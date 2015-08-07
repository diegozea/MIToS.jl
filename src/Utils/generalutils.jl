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

"""
`get_n_words(line::ASCIIString, n::Int)` returns a `Vector{ASCIIString}` with the first `n` (possibles) words/fields (delimited by space, tab or newline).
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
function get_n_words(line::ASCIIString, n::Int)
  N = length(line.data)
  coords = Array(ASCIIString, n)
  spaces = UInt8[' ', '\n', '\t']
  j = 0
  start = 0
  ended = false
  if N > 1
    block = !(line.data[1] in spaces)
    if block
      start = 1
    end
    for i in 2:N-1
      if line.data[i] in spaces
        if block
          j += 1
          if j == n
            last = line.data[N] == UInt8('\n') ? N-1 : N
            @inbounds coords[j] = ASCIIString(getindex(line.data,start:last))
            ended = true
            break
          else
            @inbounds coords[j] = ASCIIString(getindex(line.data,start:i-1))
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
      last = line.data[N] == UInt8('\n') ? N-1 : N
      @inbounds coords[j] = ASCIIString(getindex(line.data,start:last))
    end
  elseif N==1
    if line.data[1] in spaces
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
