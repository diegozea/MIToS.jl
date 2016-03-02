using Base.Intrinsics

import Base: convert, ==, !=, .==, zero, show, length,
             getindex, setindex!, rand, string

# Residues
# ========

"""
Most of the **MIToS** design is created around the `Residue` bitstype. It represents the 20 natural amino acids and a GAP value to represent insertion, deletion but also missing data: ambiguous residues and non natural amino acids.
Each residue is encoded as an integer number, this allows fast indexing operation using Residues of probability or frequency matrices.

**Residue creation and conversion**

Creation and `convert`ion of `Residue`s should be treated carefully.
`Residue` is encoded as an 8 bits type similar to `Int8`, to get faster indexing using `Int(x::Residue)`.
In this way, `Int`, `Int8` and other signed integers returns the integer value encoded by the residue.
Conversions to and from `Char`s and `Uint8` are different, to use the `Char`acter representation in IO operations.

```julia

julia> alanine = Residue('A')
A

julia> Int(alanine)
1

julia> Char(alanine)
'A'

julia> UInt8(alanine) # 0x41 == 65 == 'A'
0x41

julia> for residue in res"ARNDCQEGHILKMFPSTWYV-"
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

```
"""
bitstype 8 Residue

# Conversion from/to integers
# ---------------------------
#
# This conversions are used for indexing, comparison and operations

convert(::Type{Residue}, x::Int8) = box(Residue,unbox(Int8,x))
convert(::Type{Int8}, x::Residue) = box(Int8,unbox(Residue,x))

for typ in (Int16, Int32, Int64)
    @eval convert(::Type{$typ}, x::Residue) = convert($typ,Int8(x))
    @eval convert(::Type{Residue},  x::$typ)= convert(Residue,Int8(x))
end

Residue(x) = convert(Residue, x)

# Gaps
# ----

"""
`GAP` is the character/number representation on **MIToS** for gaps (also for non standard residues).
Lowercase characters and dots are also encoded as `GAP` in conversion from `String`s and `Char`s.
This `Residue` constant is encoded as `Residue(21)`.
"""
const GAP = Residue(21)

# Conversion from/to Char/Uint8
# -----------------------------
#
# Conversions for input and output.
# Inserts (lowercase and dots) and invalid characters will be converted into gaps ('-').

const _to_char  = Char['A','R','N','D','C','Q','E','G','H','I','L',
                       'K','M','F','P','S','T','W','Y','V','-']

convert(::Type{Char}, x::Residue) = _to_char[ Int( x ) ]

const _to_uint8 = UInt8[ UInt8(ch) for ch in  _to_char ]

convert(::Type{UInt8}, x::Residue) = _to_uint8[ Int( x ) ]

const _to_residue = fill(GAP,256)

for i in 1:length( _to_char )
    _to_residue[ Int(_to_char[i]) ] = Residue(i)
end

convert(::Type{Residue}, x::UInt8) = _to_residue[ Int(x) ];
convert(::Type{Residue}, x::Char)  = _to_residue[ Int(x) ];

# Conversion from/to ASCIIString
# ------------------------------
#
# They are useful for IO and creation of Vector{Residue}

convert(::Type{Residue}, str::ASCIIString) = Residue[ Residue(char) for char in str.data ]
convert(::Type{Residue}, str::Vector{UInt8}) = convert(Vector{Residue}, str)
convert(::Type{ASCIIString}, seq::AbstractVector{Residue}) = ASCIIString(
    UInt8[ UInt8(res) for res in seq ] )

string(seq::AbstractVector{Residue}) = ascii(seq) # "AR..." instead of the standar "[A,R,..."

# For convert a MSA stored as Vector{ASCIIString} to Matrix{Residue}
function convert(::Type{Matrix{Residue}}, sequences::Array{ASCIIString,1})
    nseq = length(sequences)
    nres = length(sequences[1])
    aln = Array(Residue,nres,nseq)
    for col in 1:nseq
        seq = sequences[col].data
        if length(seq) == nres
            aln[:,col] = Residue( seq )
        else
            throw(ErrorException(
            "There is an aligned sequence with different number of columns [ $(length(seq)) != $(nres) ]: $(ascii(seq))"
            ))
        end
    end
    aln'
end

"""
The MIToS macro `@res_str` takes a string and returns a `Vector` of `Residues` (sequence).

```julia

julia> res"MIToS"
5-element Array{MIToS.MSA.Residue,1}:
 M
 I
 T
 -
 S

```
"""
macro res_str(str)
    Residue(str)
end

# Show
# ----

show(io::IO, x::Residue) = write(io, convert(Char, x))

# Comparisons
# -----------

for fun in [:(==), :(!=), :(.==)] # , :(<), :(<=)
    @eval $(fun)(x::Residue,y::Residue) = $(fun)(Int8(x), Int8(y))
end

length(x::Residue) = length(UInt8(x))

# Random
# ------

"""
It chooses from the 20 natural residues (it doesn't generate gaps).

```julia

julia> rand(Residue)
T

julia> rand(Residue, 4, 4)
4x4 Array{MIToS.MSA.Residue,2}:
 L  E  F  L
 R  Y  L  K
 K  V  S  G
 Q  V  M  T

```
"""
rand(r, ::Type{Residue}) = Residue( rand(r, Int8(1):Int8(20)) )
