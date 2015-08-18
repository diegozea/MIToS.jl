import Base: eachline, start, done, next, eltype

type EachLineString
    string::ASCIIString
    α::Int
    len::Int
end

function eachline(string::ASCIIString)
  len = length(string)
  if len > 0
    EachLineString(string, 0, len)
  else
    throw(ErrorException("Empty String"))
  end
end

start(itr::EachLineString) = 0

done(itr::EachLineString, state) = itr.α >= itr.len || state >= itr.len

function next(itr::EachLineString, state)
  itr.α = state + 1
  index = searchindex(itr.string, '\n', state + 1)
  state = index == 0 ? itr.len : index
  itr.string[itr.α:state], state
end

eltype(::Type{EachLineString}) = ASCIIString
