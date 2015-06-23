import Base: length, getindex, size, ndims, copy, deepcopy, start, next, done, show, push!, print

# Mapping Values and Array Positions back and forth
# =================================================

type IndexedVector{T}
  values::Vector{T}
  value2index::Dict{T,Int}
end

function indexedvector{T}(vector::Vector{T})
  index = Dict{T, Int}()

  # sizehint!(index, length(vector)) # Julia 0.4
  sizehint(index, length(vector))

  i = 1
  for element in vector
    index[element] = i
    i += 1
  end

  IndexedVector(vector, index)
end

# Selection & Getters
# -------------------

selectindex{T}(indvec::IndexedVector{T}, value::T) = indvec.value2index[ value ]
selectindex{T}(indvec::IndexedVector{T}, values::Vector{T}) = Int[ indvec.value2index[value] for value in values ]
selectindex{T}(indvec::IndexedVector{T}) = selectindex(indvec, indvec.values)

selectvalue(indvec::IndexedVector) = indvec.values
selectvalue(indvec::IndexedVector, I...) = getindex(indvec.values, I...)

getindex{T,I}(indvec::IndexedVector{T}, i::AbstractArray{I} ) = indexedvector( indvec.values[ i ] )

# Copy and deepcopy
# -----------------

deepcopy(indvec::IndexedVector) = IndexedVector(deepcopy(indvec.values), deepcopy(indvec.value2index))
copy(indvec::IndexedVector) = IndexedVector(copy(indvec.values), copy(indvec.value2index))

# length, size & ndims
# --------------------

for fun in (:length, :size, :ndims)
  @eval $(fun)(indvec::IndexedVector) = $(fun)(indvec.values)
end

# Iteration
# ---------

start(indvec::IndexedVector) = 1

function next(indvec::IndexedVector, state::Int)
  value = indvec.values[state]
  ((value, state), state+1)
end

done(indvec::IndexedVector, state::Int) = done(indvec.values, state)

# Swap two vector elements
# ------------------------

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
    throw("$value is already in the IndexedVector")
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
