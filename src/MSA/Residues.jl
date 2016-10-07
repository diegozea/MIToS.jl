# Residues
# ========

if Int === Int64
    bitstype 64 Residue
else
    bitstype 32 Residue
end

@doc """
Most of the **MIToS** design is created around the `Residue` bitstype. It represents the 20
natural amino acids and a GAP value to represent insertion, deletion but also missing data:
ambiguous residues and non natural amino acids.
Each residue is encoded as an integer number, this allows fast indexing operation using
Residues of probability or frequency matrices.

**Residue creation and conversion**

Creation and `convert`ion of `Residue`s should be treated carefully.
`Residue` is encoded as an 8 bits type similar to `Int8`, to get faster indexing
using `Int(x::Residue)`. In this way, `Int`, `Int8` and other signed integers returns the
integer value encoded by the residue. Conversions to and from `Char`s and `Uint8` are
different, to use the `Char`acter representation in IO operations.

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
""" Residue

# Conversion from/to integers
# ---------------------------
#
# This conversions are used for indexing, comparison and operations

Base.convert(::Type{Residue}, x::Int) = reinterpret(Residue,x)
Base.convert(::Type{Int}, x::Residue) = reinterpret(Int,x)

# Gaps
# ----

"""
`GAP` is the character/number representation on **MIToS** for gaps (also for non standard residues).
Lowercase characters and dots are also encoded as `GAP` in conversion from `String`s and `Char`s.
This `Residue` constant is encoded as `Residue(21)`.
"""
const GAP = Residue(21) # - .

const XAA = Residue(22) # X x

# Type boundaries
# ---------------

Base.typemin(::Type{Residue}) = Residue(1)
Base.typemax(::Type{Residue}) = XAA

# Conversion from/to Char/Uint8
# -----------------------------
#
# Conversions for input and output.
# Inserts (lowercase and dots) and invalid characters will be converted into gaps ('-').

const _to_char  = Char['A','R','N','D','C','Q','E','G','H','I','L',
                       'K','M','F','P','S','T','W','Y','V','-','X' ]

function Base.convert(::Type{Char}, res::Residue)
    @inbounds char = _to_char[ Int( res ) ]
    char
end

const _max_char = 256

const _to_residue = fill(GAP, _max_char)

for index in 1:length(_to_char)
    _to_residue[ Int(_to_char[index]) ] = Residue(index)
end

for ambiguous in [  'U',   # Selenocysteine                 Sec
                    'O',   # Pyrrolysine                    Pyl
                    'B',   # Asparagine or aspartic acid    Asx     D N
                    'Z',   # Glutamine or glutamic acid     Glx     E Q
                    'J'  ] # Leucine or Isoleucine          Xle     I L

    _to_residue[ Int(ambiguous) ] = XAA
end

function Base.convert(::Type{Residue}, char::Char)
    i = Int(char)
    if i <= _max_char
        @inbounds res = _to_residue[ i ]
    else
        return GAP
    end
    res
end

# Show
# ----

function Base.show(io::IO, res::Residue)
    write(io, Char(res))
    nothing
end

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
    buffer = IOBuffer()
    for res in seq
        print(buffer, res)
    end
    takebuf_string(buffer)
end

# For convert a MSA stored as Vector{String} to Matrix{Residue}
function Base.convert(::Type{Matrix{Residue}}, sequences::Array{String,1})
    nseq = length(sequences)
    if nseq == 1
        throw(ErrorException("There are not sequences."))
    end
    nres = length(sequences[1])
    # throw() can be used with @threads : https://github.com/JuliaLang/julia/issues/17532
    for seq in sequences
        if length(seq) != nres
            throw(ErrorException(String(
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
Base.rand(r, ::Type{Residue}) = Residue(rand(r, 1:20))

# Three letters (for PDB)
# =======================

const _res2three = [ "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
                     "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
                     "   ", "XAA" ]

# Thanks to Elin- for this list
const _three2res = Dict{String, Residue}(
			                 "GLY"=>'G',
                             "ALA"=>'A',
                             "LEU"=>'L',
                             "MET"=>'M',
                             "PHE"=>'F',
                             "TRP"=>'W',
                             "LYS"=>'K',
                             "GLN"=>'Q',
                             "GLU"=>'E',
                             "SER"=>'S',
                             "PRO"=>'P',
                             "VAL"=>'V',
                             "ILE"=>'I',
                             "CYS"=>'C',
                             "TYR"=>'Y',
                             "HIS"=>'H',
                             "ARG"=>'R',
                             "ASN"=>'N',
                             "ASP"=>'D',
                             "THR"=>'T',
                             "MSE"=>'M',
                             "CSS"=>'C',
                             "2AS"=>'D',
                             "3AH"=>'H',
                             "5HP"=>'E',
                             "ACL"=>'R',
                             "AIB"=>'A',
                             "ALM"=>'A',
                             "ALO"=>'T',
                             "ALY"=>'K',
                             "ARM"=>'R',
                             "ASA"=>'D',
                             "ASB"=>'D',
                             "ASK"=>'D',
                             "ASL"=>'D',
                             "ASQ"=>'D',
                             "AYA"=>'A',
                             "BCS"=>'C',
                             "BHD"=>'D',
                             "BMT"=>'T',
                             "BNN"=>'A',
                             "BUC"=>'C',
                             "BUG"=>'L',
                             "C5C"=>'C',
                             "C6C"=>'C',
                             "CCS"=>'C',
                             "CEA"=>'C',
                             "CHG"=>'A',
                             "CLE"=>'L',
                             "CME"=>'C',
                             "CSD"=>'A',
                             "CSO"=>'C',
                             "CSP"=>'C',
                             "CSS"=>'C',
                             "CSW"=>'C',
                             "CXM"=>'M',
                             "CY1"=>'C',
                             "CY3"=>'C',
                             "CYG"=>'C',
                             "CYM"=>'C',
                             "CYQ"=>'C',
                             "DAH"=>'F',
                             "DAL"=>'A',
                             "DAR"=>'R',
                             "DAS"=>'D',
                             "DCY"=>'C',
                             "DGL"=>'E',
                             "DGN"=>'Q',
                             "DHA"=>'A',
                             "DHI"=>'H',
                             "DIL"=>'I',
                             "DIV"=>'V',
                             "DLE"=>'L',
                             "DLY"=>'K',
                             "DNP"=>'A',
                             "DPN"=>'F',
                             "DPR"=>'P',
                             "DSN"=>'S',
                             "DSP"=>'D',
                             "DTH"=>'T',
                             "DTR"=>'W',
                             "DTY"=>'Y',
                             "DVA"=>'V',
                             "EFC"=>'C',
                             "FLA"=>'A',
                             "FME"=>'M',
                             "GGL"=>'E',
                             "GLZ"=>'G',
                             "GMA"=>'E',
                             "GSC"=>'G',
                             "HAC"=>'A',
                             "HAR"=>'R',
                             "HIC"=>'H',
                             "HIP"=>'H',
                             "HMR"=>'R',
                             "HPQ"=>'F',
                             "HTR"=>'W',
                             "HYP"=>'P',
                             "IIL"=>'I',
                             "IYR"=>'Y',
                             "KCX"=>'K',
                             "LLP"=>'K',
                             "LLY"=>'K',
                             "LTR"=>'W',
                             "LYM"=>'K',
                             "LYZ"=>'K',
                             "MAA"=>'A',
                             "MEN"=>'N',
                             "MHS"=>'H',
                             "MIS"=>'S',
                             "MLE"=>'L',
                             "MPQ"=>'G',
                             "MSA"=>'G',
                             "MSE"=>'M',
                             "MVA"=>'V',
                             "NEM"=>'H',
                             "NEP"=>'H',
                             "NLE"=>'L',
                             "NLN"=>'L',
                             "NLP"=>'L',
                             "NMC"=>'G',
                             "OAS"=>'S',
                             "OCS"=>'C',
                             "OMT"=>'M',
                             "PAQ"=>'Y',
                             "PCA"=>'E',
                             "PEC"=>'C',
                             "PHI"=>'F',
                             "PHL"=>'F',
                             "PR3"=>'C',
                             "PRR"=>'A',
                             "PTR"=>'Y',
                             "SAC"=>'S',
                             "SAR"=>'G',
                             "SCH"=>'C',
                             "SCS"=>'C',
                             "SCY"=>'C',
                             "SEL"=>'S',
                             "SEP"=>'S',
                             "SET"=>'S',
                             "SHC"=>'C',
                             "SHR"=>'K',
                             "SOC"=>'C',
                             "STY"=>'Y',
                             "SVA"=>'S',
                             "TIH"=>'A',
                             "TPL"=>'W',
                             "TPO"=>'T',
                             "TPQ"=>'A',
                             "TRG"=>'K',
                             "TRO"=>'W',
                             "TYB"=>'Y',
                             "TYQ"=>'Y',
                             "TYS"=>'Y',
                             "TYY"=>'Y',
                             "AGM"=>'R',
                             "GL3"=>'G',
                             "SMC"=>'C',
                             "ASX"=>'B',
                             "CGU"=>'E',
                             "CSX"=>'C',
                             "GLX"=>'Z',
                             "UNK"=>'X',
                             "XAA"=>'X'
                            )

"""
This function returns the three letter name of the `Residue`.

```julia
julia> residue2three(Residue('G'))
"GLY"

```
"""
function residue2three(res::Residue)
    if res != GAP
        @inbounds name = _res2three[Int(res)]
    else
        throw(ErrorException("Gap has not three letter name."))
    end
    name
end

"""
It takes a three letter residue name and returns the `Residue`.

```julia
julia> three2residue("ALA")
A

```
"""
function three2residue(res::String)
    if length(res) == 3
        get(_three2res, uppercase(res), XAA)
    else
        throw(ErrorException("The residue name should have 3 letters."))
    end
end
