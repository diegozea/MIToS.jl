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
                             # Added MIToS 2.0 (XAA):
                             "UNK"=>'X',
                             "XAA"=>'X',
                             "SEC"=>'U', # Selenocysteine
                             "PYL"=>'O', # Pyrrolysine
                             "XLE"=>'J'  # Leucine or Isoleucine
                            )

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
