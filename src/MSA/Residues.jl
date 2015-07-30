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
    @eval convert(::Type{$typ}, x::Residue) = convert($typ,Int8(x))
    @eval convert(::Type{Residue},  x::$typ)= convert(Residue,Int8(x))
end

Residue(x) = convert(Residue, x)
Residue{T,N}(A::AbstractArray{T,N}) = convert(AbstractArray{Residue,N}, A)

# Gaps
# ----

const GAP = Residue(21)

# Conversion from/to Char/Uint8
# -----------------------------
#
# Conversions for input and output.
# Inserts (lowercase and dots) and invalid characters will be converted into gaps ('-').

const _to_char  = Char['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','-']

convert(::Type{Char}, x::Residue) = _to_char[ Int( x ) ]

const _to_uint8 = UInt8[ UInt8(ch) for ch in  _to_char ]

convert(::Type{UInt8}, x::Residue) = _to_uint8[ Int( x ) ]

const _to_residue = fill(GAP,256)

for i in 1:length( _to_char )
  _to_residue[ Int(_to_char[i]) ] = Residue(i)
end

convert(::Type{Residue}, x::UInt8) = _to_residue[ Int(x) ];
convert(::Type{Residue}, x::Char)  = _to_residue[ Int(x) ];

# Show
# ----

show(io::IO, x::Residue) = write(io, convert(Char, x))

# Comparisons
# -----------

for fun in [:(==),:(!=),:(<),:(<=),:(.==)]
  @eval $(fun)(x::Residue,y::Residue) = $(fun)(Int8(x), Int8(y))
end

length(x::Residue) = length(UInt8(x))
