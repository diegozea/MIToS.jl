import Base: length, getindex, size, ndims, copy, deepcopy, start, next, done, eltype, show, push!, print, ==, !=

# Mapping Values and Array Positions back and forth
# =================================================

"""
`IndexedVector{T}` allows to map values (`selectvalue`) and array positions (`selectindex`) back and forth.
Useful for vector containing unique identifiers (i.e. sequence ids in a multiple sequence alignment) that are going to be used for indexing.
Iteration over an `IndexedVector{T}` gives you `Tuple{T, Int}`.
```
julia> pdbs = IndexedVector(["2HHB", "1IRD", "1DN3"])
MIToS.MSA.IndexedVector{ASCIIString}
  values: Array(ASCIIString,(3,)) ASCIIString["2HHB","1IRD","1DN3"]
  value2index: Dict{ASCIIString,Int64} len 3
    2HHB: Int64 1
    1IRD: Int64 2
    1DN3: Int64 3

julia> for element in pdbs
          println(element)
       end
("2HHB",1)
("1IRD",2)
("1DN3",3)

julia> selectvalue(pdb_list, 2)
"1IRD"

julia> selectindex(pdb_list, "1IRD")
2
```
"""
type IndexedVector{T}
  values::Vector{T}
  value2index::Dict{T,Int}
end

"""
You can create an `IndexedVector` for any vector with unique values.
"""
function IndexedVector{T}(vector::Vector{T})
  index = Dict{T, Int}()

  sizehint!(index, length(vector))

  i = 1
  for element in vector
    if element in keys(index)
      throw(ErrorException("$element is more than one time in the vector."))
    else
      index[element] = i
      i += 1
    end
  end

  IndexedVector(vector, index)
end



# Selection & Getters
# -------------------

"""
`selectindex(indvec::IndexedVector, value)` retrieves the index of the `value` on the `IndexedVector`.
It's possible to call this method with a vector (group) of values.
This function returns always an `Int` or `Vector{Int}`.
If the value isn't in the  `IndexedVector`, `selectindex` throws a `KeyError`.
"""
selectindex{T}(indvec::IndexedVector{T}, value::T) = indvec.value2index[ value ]
selectindex{T}(indvec::IndexedVector{T}, values::Vector{T}) = Int[ indvec.value2index[value] for value in values ]
selectindex{T}(indvec::IndexedVector{T}) = Int[ i for i in 1:length(indvec)]

"""
`selectvalue(indvec::IndexedVector, index)` retrieves the value of the `index` on the `IndexedVector`.
It's possible to call this method with a vector (group) of indices or with a `Vector{Bool}` or `BitArray{1}`.
If the index isn't in the  `IndexedVector`, `selectvalue` throws a `BoundsError`.

`selectvalue{T}(indvec::IndexedVector{T}...` always return a value of type `T` or a `Vector{T}`.
You can use `getindex` (with an `Array`) in order to generate a new `IndexedVector` for the selected indices:
```
julia> pdbs = IndexedVector(["2HHB","1IRD","1DN3"]);

julia> selectvalue(pdbs, [2,3])
2-element Array{ASCIIString,1}:
 "1IRD"
 "1DN3"

julia> pdbs[[2,3]]
MIToS.MSA.IndexedVector{ASCIIString}
  values: Array(ASCIIString,(2,)) ASCIIString["1IRD","1DN3"]
  value2index: Dict{ASCIIString,Int64} len 2
    1IRD: Int64 1
    1DN3: Int64 2
```
"""
selectvalue(indvec::IndexedVector) = indvec.values
selectvalue(indvec::IndexedVector, I...) = getindex(indvec.values, I...)

getindex{T,I}(indvec::IndexedVector{T}, i::AbstractArray{I} ) = IndexedVector( indvec.values[ i ] )

# Copy and deepcopy
# -----------------

deepcopy(indvec::IndexedVector) = IndexedVector(deepcopy(indvec.values), deepcopy(indvec.value2index))
copy(indvec::IndexedVector) = IndexedVector(copy(indvec.values), copy(indvec.value2index))

# Comparisons
# -----------

for fun in (:(==), :(!==))
  @eval $(fun)(a::IndexedVector, b::IndexedVector) = $(fun)(a.values, b.values)
end

# length, size & ndims
# --------------------

for fun in (:length, :size, :ndims)
  @eval $(fun)(indvec::IndexedVector) = $(fun)(indvec.values)
end

# Iteration
# ---------

start(indvec::IndexedVector) = 1

function next{T}(indvec::IndexedVector{T}, state::Int)
  value = indvec.values[state]
  ((value, state), state+1)
end

done(indvec::IndexedVector, state::Int) = done(indvec.values, state)

eltype{T}(indvec::IndexedVector{T}) = Tuple{T,Int}

# Swap two vector elements
# ------------------------

"""
`swap!(indvec::IndexedVector, to::Int, from::Int)` interchange/swap the values on the indices `to` and `from` in the `IndexedVector`
"""
function swap!(indvec::IndexedVector, to::Int, from::Int)
  previous_id  = selectvalue(indvec, to)
  future_id    = selectvalue(indvec, from)

  indvec.values[to]   = future_id
  indvec.values[from] = previous_id

  indvec.value2index[previous_id] = from
  indvec.value2index[future_id]   = to

  indvec
end

# Push!
# =====

function push!{T}(indvec::IndexedVector{T}, value::T)
  if value in keys(indvec.value2index)
    throw(ErrorException("$value is already in the IndexedVector"))
  else
    push!(indvec.values, value)
    indvec.value2index[value] = length(indvec.values)
  end
  indvec
end

# Show & Print
# ------------

print(io::IO, indvec::IndexedVector) = dump(io, indvec)
show(io::IO, indvec::IndexedVector) = dump(io, indvec)
