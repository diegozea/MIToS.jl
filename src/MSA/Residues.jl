using Base.Intrinsics

import Base: convert, ==, !=, <, <=, .==, zero, show, length, getindex, setindex!

# Residues
# ========

bitstype 8 Residue

# Conversion from/to integers
# ---------------------------
#
# This conversions are used for indexing, comparison and operations
# They should not be used for creating Residues!

convert(::Type{Residue}, x::Int8) = box(Residue,unbox(Int8,x))
convert(::Type{Int8}, x::Residue) = box(Int8,unbox(Residue,x))

for typ in (Int16, Int32, Int64)
    @eval convert(::Type{$typ}, x::Residue) = convert($typ,int8(x))
    @eval convert(::Type{Residue},  x::$typ)= convert(Residue,int8(x))
end

residue(x) = convert(Residue, x)
residue{T,N}(A::AbstractArray{T,N}) = convert(AbstractArray{Residue,N}, A)

# Gaps
# ----
#
# Gap is 0x00 (zero). This is useful for using Sparse matrix or vector.

const GAP = residue(0)

zero(::Type{Residue}) = GAP

# Conversion from/to Char/Uint8
# -----------------------------
#
# Conversions for input and output.
# Inserts (lowercase and dots) and invalid characters will be converted into gaps ('-').

const _to_char  = Char['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']

convert(::Type{Char}, x::Residue) = x == GAP ? '-' : _to_char[ int( x ) ]

const _to_uint8 = uint8( _to_char )

convert(::Type{Uint8}, x::Residue) = x == GAP ? 0x2d : _to_uint8[ int( x ) ]

const _to_residue = zeros(Residue,256)

for i in 1:length( _to_char )
  _to_residue[ _to_char[i] ] = residue(i)
end

convert(::Type{Residue}, x::Uint8) = _to_residue[ int(x) ];
convert(::Type{Residue}, x::Char)  = _to_residue[ int(x) ];

# Show
# ----

show(io::IO, x::Residue) = write(io, convert(Char, x))

# Comparisons
# -----------

for fun in [:(==),:(!=),:(<),:(<=),:(.==)]
  @eval $(fun)(x::Residue,y::Residue) = $(fun)(int8(x), int8(y))
end

length(x::Residue) = length(uint8(x))
