# Three letters (for PDB)
# =======================

const _res2three = [ "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
                     "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
                     "   ", "XAA" ]

# Thanks to Elin- for this list
const _three2res = Dict{String, Residue}(
    "BMT"=>'T',
    "DCY"=>'C',
    "SCY"=>'C',
    "ASB"=>'D',
    "STY"=>'Y',
    "GLY"=>'G',
    "ASN"=>'N',
    "ASA"=>'D',
    "DTR"=>'W',
    "ARM"=>'R',
    "FLA"=>'A',
    "CYQ"=>'C',
    "DSP"=>'D',
    "SHC"=>'C',
    "BHD"=>'D',
    "DVA"=>'V',
    "CSD"=>'A',
    "PHI"=>'F',
    "PHE"=>'F',
    "TYB"=>'Y',
    "HIP"=>'H',
    "DIL"=>'I',
    "HYP"=>'P',
    "TYR"=>'Y',
    "IYR"=>'Y',
    "3AH"=>'H',
    "NLE"=>'L',
    "SOC"=>'C',
    "SVA"=>'S',
    "CSW"=>'C',
    "LYZ"=>'K',
    "TPL"=>'W',
    "ALY"=>'K',
    "GLU"=>'E',
    "MSE"=>'M',
    "PCA"=>'E',
    "DHI"=>'H',
    "5HP"=>'E',
    "C5C"=>'C',
    "GLN"=>'Q',
    "CY3"=>'C',
    "TYS"=>'Y',
    "AGM"=>'R',
    "LYS"=>'K',
    "MEN"=>'N',
    "DPR"=>'P',
    "TRP"=>'W',
    "THR"=>'T',
    "TYQ"=>'Y',
    "VAL"=>'V',
    "CEA"=>'C',
    "DAL"=>'A',
    "MSA"=>'G',
    "DLY"=>'K',
    "BUC"=>'C',
    "OAS"=>'S',
    "NMC"=>'G',
    "TPQ"=>'A',
    "DLE"=>'L',
    "BNN"=>'A',
    "AYA"=>'A',
    "DNP"=>'A',
    "PAQ"=>'Y',
    "DTY"=>'Y',
    "CCS"=>'C',
    "AIB"=>'A',
    "SAR"=>'G',
    "ARG"=>'R',
    "PTR"=>'Y',
    "SHR"=>'K',
    "FME"=>'M',
    "MIS"=>'S',
    "NLN"=>'L',
    "NLP"=>'L',
    "PHL"=>'F',
    "TYY"=>'Y',
    "CSP"=>'C',
    "CSS"=>'C',
    "MET"=>'M',
    "BCS"=>'C',
    "MAA"=>'A',
    "PR3"=>'C',
    "TPO"=>'T',
    "GLZ"=>'G',
    "HMR"=>'R',
    "ASK"=>'D',
    "DAH"=>'F',
    "MLE"=>'L',
    "NEM"=>'H',
    "CHG"=>'A',
    "DGL"=>'E',
    "DTH"=>'T',
    "OCS"=>'C',
    "GMA"=>'E',
    "TRO"=>'W',
    "CLE"=>'L',
    "ALO"=>'T',
    "GL3"=>'G',
    "BUG"=>'L',
    "CXM"=>'M',
    "HPQ"=>'F',
    "LLY"=>'K',
    "SCH"=>'C',
    "DHA"=>'A',
    "TRG"=>'K',
    "CME"=>'C',
    "SER"=>'S',
    "ALM"=>'A',
    "MPQ"=>'G',
    "SMC"=>'C',
    "CSO"=>'C',
    "LTR"=>'W',
    "ASP"=>'D',
    "DAS"=>'D',
    "GSC"=>'G',
    "OMT"=>'M',
    "SCS"=>'C',
    "HIC"=>'H',
    "DSN"=>'S',
    "DGN"=>'Q',
    "TIH"=>'A',
    "ALA"=>'A',
    "KCX"=>'K',
    "IIL"=>'I',
    "NEP"=>'H',
    "CYS"=>'C',
    "GGL"=>'E',
    "LLP"=>'K',
    "HAC"=>'A',
    "DIV"=>'V',
    "C6C"=>'C',
    "DPN"=>'F',
    "SAC"=>'S',
    "CY1"=>'C',
    "CYM"=>'C',
    "SEL"=>'S',
    "EFC"=>'C',
    "DAR"=>'R',
    "SET"=>'S',
    "ASL"=>'D',
    "HIS"=>'H',
    "SEP"=>'S',
    "MVA"=>'V',
    "PRO"=>'P',
    "ASQ"=>'D',
    "PRR"=>'A',
    "CGU"=>'E',
    "CSX"=>'C',
    "MHS"=>'H',
    "ILE"=>'I',
    "ACL"=>'R',
    "HTR"=>'W',
    "PEC"=>'C',
    "LEU"=>'L',
    "LYM"=>'K',
    "2AS"=>'D',
    "CYG"=>'C',
    "HAR"=>'R' )

"""
This function returns the three letter name of the `Residue`.

```julia
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
