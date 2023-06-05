const _res2three = [ "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
                     "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
                     "   ", "XAA" ]

"""
This function returns the three letter name of the `Residue`.

```jldoctest
julia> using MIToS.MSA

julia> residue2three(Residue('G'))
"GLY"

```
"""
function residue2three(res::Residue)
    int_res = Int(res)
    if int_res == 21  || !isvalid(res)
        #         GAP
        throw(ErrorException("Residue($(int_res)) has not three letter name."))
    end
    _res2three[int_res]
end

"""
It takes a three letter residue name and returns the corresponding `Residue`.
If the name isn't in the MIToS dictionary, a `XAA` is returned.

```jldoctest
julia> using MIToS.MSA

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

const _three2res = convert(Dict{String, Residue}, THREE2ONE) # THREE2ONE from MIToS.Utils
