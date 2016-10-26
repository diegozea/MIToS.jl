abstract ResidueAlphabet

immutable GappedAlphabet <: ResidueAlphabet end

immutable UngappedAlphabet <: ResidueAlphabet end

immutable ReducedAlphabet <: ResidueAlphabet
    mapping::Vector{Int}
    len::Int
end

# Creation
# --------

function (::Type{ReducedAlphabet})(str::AbstractString)
    mapping = fill!(Array(Int, Int(XAA)), Int(XAA))
    ingroup = false
    pos = 0
    for char in str
        if char == '('
            ingroup = true
            pos += 1
            continue
        end
        if char == ')'
            ingroup = false
            continue
        end
        if !ingroup
            pos += 1
        end
        mapping[Int(Residue(char))] = pos
    end
    ReducedAlphabet(mapping, pos)
end

macro reduced_str(str)
    ReducedAlphabet(str)
end

# Iteration Interface
# -------------------

Base.length(::UngappedAlphabet)  = 20
Base.length(::GappedAlphabet)    = 21
Base.length(ab::ReducedAlphabet) = ab.len

Base.start(::ResidueAlphabet) = 1

Base.next(::ResidueAlphabet, state) = (state, state + 1)

Base.done(ab::UngappedAlphabet, state) = state > length(ab)
Base.done(ab::GappedAlphabet, state)   = state > length(ab)
Base.done(ab::ReducedAlphabet, state)  = state > length(ab)

Base.eltype(::ResidueAlphabet) = Int

# Show
# ----

function Base.show(io::IO, ab::ResidueAlphabet)
    print(io, typeof(ab), " of length ", length(ab), ". Residues : res\"",
          join(_to_char[1:length(ab)]), "\"")
end

function Base.show(io::IO, ab::ReducedAlphabet)
    groups = ab.mapping
    print(io, "ReducedAlphabet of length ", length(ab), " : reduced\"")
    for i in ab
        chars = _to_char[groups .== i]
        if length(chars) == 1
            print(io, chars[1])
        elseif length(chars) > 1
            print(io, "(", join(chars),")")
        end
    end
    println(io, "\"")
end

# getindex
# --------

@inline function Base.getindex(ab::ResidueAlphabet, res::Residue)::Int
    ifelse(Int(res) <= length(ab), Int(res), 22)
end

@inline function Base.getindex(ab::ReducedAlphabet, res::Residue)::Int
    @inbounds value = ab.mapping[Int(res)]
    value
end

# In Alphabet
# -----------

Base.in(res::Residue, alphabet::ResidueAlphabet) = alphabet[res] != 22
