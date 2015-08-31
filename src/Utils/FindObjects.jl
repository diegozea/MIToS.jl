abstract AbstractTest

abstract TestType <: AbstractTest
abstract TestOperation <: AbstractTest

immutable Is{T} <: TestType
 field::Symbol
 value::T
end

immutable In{T} <: TestType
 field::Symbol
 value::T
end

_test{T}(object, test::In{T}) = getfield(object, test.field) in test.value

_test{T}(object, test::Is{T}) = getfield(object, test.field) == test.value
_test(object, test::Is{Regex}) = ismatch(test.value, getfield(object, test.field))
_test(object, test::Is{Function}) = test.value(getfield(object, test.field))

immutable Not{T <: TestType} <: TestOperation
 test::T
end

_test{T}(object, test::Not{T}) = !_test(object, test.test::T)

"Returns `true` if all the tests are satisfied"
function isobject(object, tests::AbstractTest...)
  for test in tests
    if !_test(object, test)
      return(false)
    end
  end
  true
end

"""Returns a `Nullable` of the `field`.
The object is `null` if any test fails"""
function capture(object, field::Symbol, tests::AbstractTest...)
  captured = getfield(object, field)
  for test in tests
    if !_test(object, test)
      return(Nullable{typeof(captured)}())
    end
  end
  Nullable(captured)
end

"Returns a vector of the indexes for which `isobject` is true."
function findobjects(vector, tests::AbstractTest...)
  len = length(vector)
  indexes = Array(Int, len)
  n = 0
  for i in 1:len
    if isobject(vector[i], tests...)
      n += 1
      indexes[n] = i
    end
  end
  resize!(indexes, n)
end

"Returns a vector with the objects for which `isobject` is true given the Tests."
function collectobjects{T}(vector::AbstractVector{T}, tests::AbstractTest...)
  len = length(vector)
  elements = Array(T, len)
  n = 0
  for i in 1:len
    if isobject(vector[i], tests...)
      n += 1
      elements[n] = vector[i]
    end
  end
  resize!(elements, n)
end

"Returns the **type** of the `field` in the first object of the collection"
function guess_type(collection, field::Symbol)
  for object in collection
    if field in fieldnames(object)
      return(fieldtype(typeof(object), field))
    end
  end
  Any
end

"""Returns a vector of Nullables with the captures of the `field`s.
Null if any test fails or the object hasn't the `field`."""
function collectcaptures(vector, field::Symbol, tests::AbstractTest...)
  Nullable{guess_type(vector, field)}[ field in fieldnames(object) ? capture(object, field, tests...) : Nullable() for object in vector ]
end
