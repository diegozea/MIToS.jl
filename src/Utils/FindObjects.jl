abstract AbstractTest

immutable Not{T} <: AbstractTest
 field::Symbol
 value::T
end

immutable Is{T} <: AbstractTest
 field::Symbol
 value::T
end

_test_not{T}(object, test::Is{T}) =         getfield(object, test.field) != test.value
_test_not(object,    test::Is{Regex}) =    !ismatch(test.value, getfield(object, test.field))
_test_not(object,    test::Is{Function}) = !test.value(getfield(object, test.field))

_test_not{T}(object, test::Not{T}) =        getfield(object, test.field) == test.value
_test_not(object,    test::Not{Regex}) =    ismatch(test.value, getfield(object, test.field))
_test_not(object,    test::Not{Function}) = test.value(getfield(object, test.field))

function isobject(object, tests...)
  for test in tests
    if _test_not(object, test)
      return(false)
    end
  end
  true
end

"Returns a vector of the indexes where `isobject` is true."
function findobjects(vector, tests...)
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