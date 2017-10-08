"Abstract type to define residue alphabet types."
abstract type ResidueAlphabet end

"""
This type defines the usual alphabet of the 20 natural residues and a gap character.

```julia
julia> GappedAlphabet()
MIToS.MSA.GappedAlphabet of length 21. Residues : res"ARNDCQEGHILKMFPSTWYV-"

```
"""
struct GappedAlphabet <: ResidueAlphabet end

"""
This type defines the usual alphabet of the 20 natural residues, without the gap character.

```julia
julia> UngappedAlphabet()
MIToS.MSA.UngappedAlphabet of length 20. Residues : res"ARNDCQEGHILKMFPSTWYV"

```
"""
struct UngappedAlphabet <: ResidueAlphabet end

"""
`ReducedAlphabet` allows the construction of reduced residue alphabets, where residues
inside parenthesis belong to the same group.

```julia
julia> ab = ReducedAlphabet("(AILMV)(RHK)(NQST)(DE)(FWY)CGP")
MIToS.MSA.ReducedAlphabet of length 8 : "(AILMV)(RHK)(NQST)(DE)(FWY)CGP"

julia> ab[Residue('K')]
2

```
"""
struct ReducedAlphabet <: ResidueAlphabet
    mapping::Array{Int,1} # Residue -> Int
    named::NamedArray{Int,1,Array{Int,1},Tuple{OrderedDict{String,Int}}}
    len::Int
end

# Creation
# --------

function ReducedAlphabet(str::AbstractString)
    N = Int(XAA)
    mapping = fill!(Array(Int, N), N)
    ingroup = false
    pos = 0
    group_names = fill!(Array(String, N), "")
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
        int_residue = Int(Residue(char))
        @assert int_residue != 22 "$char isn't valid for a residue alphabet." # N == 22
        group_names[pos] = string(group_names[pos], char)
        mapping[int_residue] = pos
    end
    named = NamedArray(collect(1:pos), (
        OrderedDict{String,Int}(group_names[i] => i for i in 1:pos),), ("Groups",))
    ReducedAlphabet(mapping, named, pos)
end

# Iteration Interface
# -------------------

Base.length(ab::UngappedAlphabet)  = 20
Base.length(ab::GappedAlphabet)    = 21
Base.length(ab::ReducedAlphabet) = ab.len

Base.start(ab::ResidueAlphabet) = 1

Base.next(ab::ResidueAlphabet, state) = (state, state + 1)

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
    groups = names(ab.named, 1)
    print(io, "MIToS.MSA.ReducedAlphabet of length ", length(ab), " : \"")
    for i in ab
        chars = groups[i]
        if length(chars) == 1
            print(io, chars[1])
        elseif length(chars) > 1
            print(io, "(", join(chars),")")
        end
    end
    print(io, "\"")
end

# getindex
# --------

@inline function Base.getindex(ab::ResidueAlphabet, res::Int)::Int
    ifelse(res <= length(ab), res, 22)
end

@inline function Base.getindex(ab::ReducedAlphabet, res::Int)::Int
    @inbounds value = ab.mapping[res]
    value
end

@inline function Base.getindex(ab::ResidueAlphabet, res::Residue)::Int
    ab[Int(res)]
end

@inline function Base.getindex(ab::ResidueAlphabet, res::String)::Int
    @assert length(res) == 1 "The string with the residue should have only one character."
    i = Int(Residue(res[1]))
    ifelse(i <= length(ab), i, 22)
end

@inline function Base.getindex(ab::ReducedAlphabet, res::String)::Int
    @inbounds value = ab.named[res]
    value
end


# In Alphabet
# -----------

Base.in(res::Residue, alphabet::ResidueAlphabet) = Int(res) <= length(alphabet)
Base.in(res::Residue, alphabet::ReducedAlphabet) = alphabet[res] != 22

# Names
# -----

"""
It returns the name of each group. The name is a string with the one letter code of each
residue that belong to the group.

```julia
julia> ab = ReducedAlphabet("(AILMV)(RHK)(NQST)(DE)(FWY)CGP")
MIToS.MSA.ReducedAlphabet of length 8 : "(AILMV)(RHK)(NQST)(DE)(FWY)CGP"

julia> names(ab)
8-element Array{String,1}:
 "AILMV"
 "RHK"
 "NQST"
 "DE"
 "FWY"
 "C"    
 "G"
 "P"

```
"""
Base.names(alphabet::ReducedAlphabet) = names(alphabet.named, 1)
Base.names(alphabet::ResidueAlphabet) = String[ string(Residue(i)) for i in alphabet ]

# Dict of names to indexes
# ------------------------

@inline _getdict(n::NamedArray{Int,1,Array{Int,1},Tuple{OrderedDict{String,Int}}}) = n.dicts[1]

const _UngappedAlphabet_Names = OrderedDict{String,Int}(string(Residue(i))=>i for i in 1:20)
const _GappedAlphabet_Names = OrderedDict{String,Int}(string(Residue(i))=>i for i in 1:21)

@inline getnamedict(alphabet::UngappedAlphabet) = _UngappedAlphabet_Names
@inline getnamedict(alphabet::GappedAlphabet) = _GappedAlphabet_Names

"""
It takes a `ResidueAlphabet` and returns a dictionary from group name to group position.

```julia
julia> ab = ReducedAlphabet("(AILMV)(RHK)(NQST)(DE)(FWY)CGP")
MIToS.MSA.ReducedAlphabet of length 8 : "(AILMV)(RHK)(NQST)(DE)(FWY)CGP"

julia> getnamedict(ab)
DataStructures.OrderedDict{String,Int64} with 8 entries:
  "AILMV" => 1
  "RHK"   => 2
  "NQST"  => 3
  "DE"    => 4
  "FWY"   => 5
  "C"     => 6
  "G"     => 7
  "P"     => 8

```
"""
@inline getnamedict(alphabet::ReducedAlphabet) = _getdict(copy(alphabet.named))
