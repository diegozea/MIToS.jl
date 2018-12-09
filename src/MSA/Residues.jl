# Residues
# ========

if Int === Int64
    primitive type Residue 64 end
else
    primitive type Residue 32 end
end

@doc """
Most of the **MIToS** design is created around the `Residue` bitstype. It has
representations for the 20 natural amino acids, a value representing insertions and
deletions (`GAP`, `'-'`) and one representing unknown, ambiguous and non standard residues
(`XAA`, `'X'`).
Each `Residue` is encoded as an integer number, with the same bit representation and size
than a `Int`. This allows fast indexing operation of probability or frequency matrices.

**Residue creation and conversion**

Creation and conversion of `Residue`s should be treated carefully. `Residue` is encoded
as a 32 or 64 bits type similar to `Int`, to get fast indexing using `Int(x::Residue)`.
`Int` simply calls `reinterpret` without checking if the residue is valid.
Valid residues have integer values in the closed interval [1,22]. `convert` from `Int`  and
`Char` always returns valid residues, however it's possible to find invalid residues
(they are shown using the character `'�'`) after the creation of uninitialized `Residue`
arrays (i.e. using `Array`). You can use `zeros`, `ones` or `rand` to get initialized
`Residue` arrays with valid residues.
Conversions to and from `Char`s changes the bit representation and allows the use of the
usual character representation of residues and amino acids. This conversions are used in IO
operations and always return valid residues. In conversions from `Char`, lowercase letters,
`'*'`, `'-'` and `'.'` are translated to `GAP`, letters representing the 20 natural amino
(ARNDCQEGHILKMFPSTWYV) acids are translated to their corresponding `Residue` and any other
character is translated to `XAA`. Since lowercase letters and dots are translated to gaps,
Pfam MSA insert columns are converted to columns full of gaps.

```jldoctest
julia> using MIToS.MSA

julia> alanine = Residue('A')
A

julia> Char(alanine)
'A': ASCII/Unicode U+0041 (category Lu: Letter, uppercase)

julia> for residue in res"ARNDCQEGHILKMFPSTWYV-X"
           println(residue, " ", Int(residue))
       end
A 1
R 2
N 3
D 4
C 5
Q 6
E 7
G 8
H 9
I 10
L 11
K 12
M 13
F 14
P 15
S 16
T 17
W 18
Y 19
V 20
- 21
X 22

```
""" Residue

# Conversion from/to integers
# ---------------------------
#
# This conversions are used for indexing, comparison and operations

"""
It takes an `Int` and returns the `Int` value used to represent a valid `Residue`.
Invalid residues are encoded as the integer 23.
"""
@inline _valid_residue_integer(x::Int) = ifelse(1 <= x <= 22, x, 22)
#                                                                XAA

Base.convert(::Type{Residue}, x::Int) = reinterpret(Residue, _valid_residue_integer(x))

# Conversion to `Int` doesn’t check if the residue is valid
@inline Base.convert(::Type{Int}, x::Residue) = reinterpret(Int,x)

Residue(i::Int) = convert(Residue, i)
Base.Int(res::Residue) = convert(Int, res)

# ndims
# -----

Base.ndims(r::Residue) = 0
Base.ndims(::Type{Residue}) = 0

# Scalar
# ------

Base.broadcastable(x::Residue) = Ref(x)

# Gaps
# ----

"""
`GAP` is the `Residue` representation on MIToS for gaps (`'-'`, insertions and
deletions). Lowercase residue characters, dots and `'*'` are encoded as `GAP` in conversion
from `String`s and `Char`s. This `Residue` constant is encoded as `Residue(21)`.
"""
const GAP = Residue(21) # - .

"""
`XAA` is the `Residue` representation for unknown, ambiguous and non standard residues.
This `Residue` constant is encoded as `Residue(22)`.
"""
const XAA = Residue(22) # X

# Check for valid residues
# ------------------------

"""
`isvalid(res::Residue)`

It returns `true` if the encoded integer is in the closed interval [1,22].
"""
Base.isvalid(::Type{Residue}, res::Residue) = 1 <= reinterpret(Int,res) <= 22
Base.isvalid(res::Residue) = isvalid(Residue, res)

# Type boundaries
# ---------------

Base.typemin(::Type{Residue}) = Residue(1)
Base.typemax(::Type{Residue}) = Residue(22)

# zeros and ones
# --------------

Base.zero(::Type{Residue}) = GAP
Base.one( ::Type{Residue}) = XAA

# Conversion from/to Char/Uint8
# -----------------------------
#
# Conversions for input and output.

const _to_char  = Char['A','R','N','D','C','Q','E','G','H','I','L',
                       'K','M','F','P','S','T','W','Y','V','-','X']

function Base.convert(::Type{Char}, res::Residue)
    @inbounds char = _to_char[ _valid_residue_integer(Int(res)) ]
    char
end

"""
'z' is the maximum between 'A':'Z', 'a':'z', '.', '-' and '*'.
'z' is 'GAP' but the next character to 'z' is '{', i.e. `XAA`.
"""
const _max_char = Int('z') + 1

const _to_residue = fill(XAA, _max_char)

# Residues
for index in 1:22
    _to_residue[ Int(_to_char[index]) ] = Residue(index)
end

# A lowercase Char is translated to GAP
for char in 'a':'z'
    _to_residue[ Int(char) ] = GAP
end

_to_residue[ Int('.') ] = GAP # Gap in insert columns
_to_residue[ Int('-') ] = GAP # Usual GAP character
_to_residue[ Int('*') ] = GAP # Usual representation of a translated stop codon

@inline Base.@pure function Base.convert(::Type{Residue}, char::Char)::Residue
    @inbounds if char < '{'
        _to_residue[ Int(char) ]
    else
        XAA
    end
end

Base.Char(res::Residue) = convert(Char, res)
Residue(char::Char) = convert(Residue, char)

# Bits
# ----

Base.bitstring(res::Residue) = bitstring(reinterpret(Int, res))

# Show
# ----

# Invalid residues are shown with the '�' character
function _write(io::IO, res::Residue)
    if isvalid(res)
        @inbounds char = _to_char[Int(res)]
        write(io, char)
    else
        write(io, '�')
    end
    nothing
end

Base.show(io::IO, res::Residue) = _write(io, res)

# Print
# -----

# Invalid residues are printed with the 'X' character
function Base.print(io::IO, res::Residue)
    write(io, Char(res))
    nothing
end

# Conversion from/to String
# ------------------------------
#
# They are useful for IO and creation of Vector{Residue}

Base.convert(::Type{Vector{Residue}}, str::AbstractString) = Residue[ char for char in str ]

function Base.convert(::Type{String}, seq::Vector{Residue})
    # Buffer length can be length(seq) since Char(res) is always ASCII
    buffer = IOBuffer(Array{UInt8}(undef, length(seq)), read=true, write=true)
    # To start at the beginning of the buffer:
    truncate(buffer,0)
    for res in seq
        write(buffer, Char(res))
    end
    String(take!(buffer))
end

Base.String(seq::Vector{Residue}) = convert(String, seq)

function _get_msa_size(sequences::Array{String,1})
    nseq = length(sequences)
    if nseq == 0
        throw(ErrorException("There are not sequences."))
    end
    nres = length(sequences[1])
    nseq, nres
end

function _convert_to_matrix_residues(sequences::Array{String,1}, size::Tuple{Int,Int})
   nseq, nres = size
   aln = Array{Residue}(undef, nseq, nres)
   @inbounds for (i, str) in enumerate(sequences)
       for (j, char) in enumerate(str)
           aln[CartesianIndex(i, j)] = Residue(char)
       end
   end
   aln
end

# For convert a MSA stored as Vector{String} to Matrix{Residue}
# This checks that all the sequences have the same length
function Base.convert(::Type{Matrix{Residue}}, sequences::Array{String,1})
    nseq, nres = _get_msa_size(sequences)
    # throw() can be used with @threads : https://github.com/JuliaLang/julia/issues/17532
    for seq in sequences
        if length(seq) != nres
            throw(ErrorException(string(
                "There is an aligned sequence with different number of columns",
                "[ ", length(seq), " != ", nres, " ]: ", String(seq) )))
        end
    end
    _convert_to_matrix_residues(sequences, (nseq, nres))
end

# Non-Standard String Literal
# ---------------------------

"""
The MIToS macro `@res_str` takes a string and returns a `Vector` of `Residues` (sequence).

```jldoctest
julia> using MIToS.MSA

julia> res"MIToS"
5-element Array{Residue,1}:
 M
 I
 T
 -
 S

```
"""
macro res_str(str)
    convert(Vector{Residue}, str)
end

# Comparisons
# -----------

Base.:(==)(x::Residue, y::Residue) = x === y
Base.:(!=)(x::Residue, y::Residue) = x !== y

Base.length(res::Residue) = length(Int(res))

# Random
# ------

"""
It chooses from the 20 natural residues (it doesn't generate gaps).

```jldoctest
julia> using MIToS.MSA

julia> using Random

julia> Random.seed!(1); # Reseed the random number generator.

julia> rand(Residue)
P

julia> rand(Residue, 4, 4)
4×4 Array{Residue,2}:
 N  N  T  D
 G  Y  L  I
 R  V  F  L
 P  C  K  L

```
"""
Random.rand(rng::AbstractRNG, ::Random.SamplerType{Residue}) = Residue(rand(rng, 1:20))

function Random.Sampler(RNG::Type{<:AbstractRNG}, res::Residue, r::Random.Repetition)
    Random.SamplerSimple(res, Random.Sampler(RNG, 1:20, r))
end

Random.rand(rng::AbstractRNG, sp::Random.SamplerSimple{Residue}) = rand(rng, sp.data)
