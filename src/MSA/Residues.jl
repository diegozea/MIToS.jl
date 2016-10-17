# Residues
# ========

if Int === Int64
    bitstype 64 Residue
else
    bitstype 32 Residue
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

```julia
julia> alanine = Residue('A')
A

julia> Char(alanine)
'A'

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

const _max_char = Int('z') # 'z' is the maximum between 'A':'Z', 'a':'z', '.', '-' and '*'

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

function Base.convert(::Type{Residue}, char::Char)
    i = Int(char)
    if 1 <= i <= _max_char
        @inbounds res = _to_residue[ i ]
    else
        res = XAA
    end
    res
end

# Bits
# ----

Base.bits(res::Residue) = bits(reinterpret(Int, res))

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
    #                 data                         readable    writable
    buffer = IOBuffer(Array(UInt8, length(seq)),    true,       true)
    # To start at the beginning of the buffer:
    truncate(buffer,0)
    for res in seq
        write(buffer, Char(res))
    end
    takebuf_string(buffer)
end

# For convert a MSA stored as Vector{String} to Matrix{Residue}
function Base.convert(::Type{Matrix{Residue}}, sequences::Array{String,1})
    nseq = length(sequences)
    if nseq == 0
        throw(ErrorException("There are not sequences."))
    end
    nres = length(sequences[1])
    # throw() can be used with @threads : https://github.com/JuliaLang/julia/issues/17532
    for seq in sequences
        if length(seq) != nres
            throw(ErrorException(string(
                "There is an aligned sequence with different number of columns",
                "[ ", length(seq), " != ", nres, " ]: ", String(seq) )))
        end
    end
    aln = Array(Residue, nseq, nres)
    # @inbounds @threads for i in 1:nseq
    @inbounds for i in 1:nseq
        aln[i,:] = collect(sequences[i])
    end
    aln
end

# Non-Standard String Literal
# ---------------------------

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
    convert(Vector{Residue}, str)
end

# Comparisons
# -----------

Base.:(==)(x::Residue, y::Residue) = Int(x) == Int(y)
Base.:(!=)(x::Residue, y::Residue) = Int(x) != Int(y)

Base.length(res::Residue) = length(Int(res))

# Random
# ------

"""
It chooses from the 20 natural residues (it doesn't generate gaps).

```julia
julia> srand(1); # Reseed the random number generator.

julia> rand(Residue)
Y

julia> rand(Residue, 4, 4)
4×4 Array{MIToS.MSA.Residue,2}:
 E  L  K  P
 V  N  A  M
 I  D  A  F
 F  Y  C  P

```
"""
Base.rand(r, ::Type{Residue}) = Residue(rand(r, 1:20))
