"""
`PDBFile <: Format`

Protein Data Bank (PDB) format.
It provides a standard representation for macromolecular
structure data derived from X-ray diffraction and NMR studies.
"""
immutable PDBFile <: Format end

"""
`parse(io, ::Type{PDBFile}; chain="all", model="all", group="all", atomname="all", onlyheavy=false, occupancyfilter=false)`

Reads a text file of a PDB entry.
Returns a list of `PDBResidue` (view `MIToS.PDB.PDBResidues`).
Setting `chain`, `model`, `group`, `atomname` and `onlyheavy` values
can be used to select of a subset of all residues. Group can be `ATOM`
or `HETATM`. If not set, all residues are returned.
If the keyword argument `occupancyfilter` (default: `false`) is `true`,
only the atoms with the best occupancy are returned.
"""
function parse(io::Union{IO, ASCIIString}, ::Type{PDBFile}; chain::ASCIIString = "all",
               model::ASCIIString = "all", group::ASCIIString = "all", atomname::ASCIIString="all",
               onlyheavy::Bool=false, occupancyfilter::Bool=false)
    residue_dict = OrderedDict{PDBResidueIdentifier, Vector{PDBAtom}}()
    atom_model = 0
    for line in eachline(io)
        line_id = length(line) < 6 ? replace(line, ' ', "") : replace(line[1:6], ' ', "") # i.e. "END\n"
        if line_id == "MODEL"
            atom_model += 1
        end
        if (group=="all" && (line_id=="ATOM" || line_id=="HETATM")) || (line_id == group)
            atom_chain = string(line[22])
            atom_name = replace(line[13:16],' ',"")
            element = replace(line[77:78],' ',"")
            if  (chain=="all" || chain==atom_chain) && ((model==0 || model=="all") || model==string(atom_model+1)) &&
                    (atomname=="all" || atomname==atom_name) && (!onlyheavy || element!="H")

                # 23 - 26        Integer         Residue sequence number.
                # 27             AChar           Code for insertion of residues.
                PDB_number = replace(line[23:27],' ',"")

                name = replace(line[18:20],' ',"")
                x = float(replace(line[31:38],' ',""))
                y = float(replace(line[39:46],' ',""))
                z = float(replace(line[47:54],' ',""))
                occupancy = float(replace(line[55:60],' ',""))
                B = replace(line[61:66],' ',"")

                mdl = atom_model == 0 ? "1" : string(atom_model)

                residue_id = PDBResidueIdentifier("", PDB_number, name, line_id, mdl, atom_chain)
                atom_data  = PDBAtom(Coordinates(x,y,z), atom_name, element, occupancy, B)

                value = get!(residue_dict, residue_id, PDBAtom[])
                push!(value, atom_data)

            end
        end
    end
    _generate_residues(residue_dict, occupancyfilter)
end

# Print PDB
# =========

# ATOM & HETATM
# COLUMNS        DATA TYPE       CONTENTS
# --------------------------------------------------------------------------------
#  1 -  6        Record name     "ATOM  "
#  7 - 11        Integer         Atom serial number.
# 13 - 16        Atom            Atom name.
# 17             Character       Alternate location indicator.
# 18 - 20        Residue name    Residue name.
# 22             Character       Chain identifier.
# 23 - 26        Integer         Residue sequence number.
# 27             AChar           Code for insertion of residues.
# 31 - 38        Real(8.3)       Orthogonal coordinates for X in Angstroms.
# 39 - 46        Real(8.3)       Orthogonal coordinates for Y in Angstroms.
# 47 - 54        Real(8.3)       Orthogonal coordinates for Z in Angstroms.
# 55 - 60        Real(6.2)       Occupancy.
# 61 - 66        Real(6.2)       Temperature factor (Default = 0.0).
# 73 - 76        LString(4)      Segment identifier, left-justified.
# 77 - 78        LString(2)      Element symbol, right-justified.
# 79 - 80        LString(2)      Charge on the atom.
# Example:
#          1         2         3         4         5         6         7         8
# 12345678901234567890123456789012345678901234567890123456789012345678901234567890
# ATOM    145  N   VAL A  25      32.433  16.336  57.540  1.00 11.92      A1   N
# ATOM    146  CA  VAL A  25      31.132  16.439  58.160  1.00 11.85      A1   C
# ATOM    147  C   VAL A  25      30.447  15.105  58.363  1.00 12.34      A1   C
# ATOM    148  O   VAL A  25      29.520  15.059  59.174  1.00 15.65      A1   O
# ATOM    149  CB AVAL A  25      30.385  17.437  57.230  0.28 13.88      A1   C
# ATOM    150  CB BVAL A  25      30.166  17.399  57.373  0.72 15.41      A1   C
# ATOM    151  CG1AVAL A  25      28.870  17.401  57.336  0.28 12.64      A1   C
# ATOM    152  CG1BVAL A  25      30.805  18.788  57.449  0.72 15.11      A1   C
# ATOM    153  CG2AVAL A  25      30.835  18.826  57.661  0.28 13.58      A1   C
# ATOM    154  CG2BVAL A  25      29.909  16.996  55.922  0.72 13.25      A1   C
#
# HETATM 1357 MG    MG   168       4.669  34.118  19.123  1.00  3.16          MG2+
# HETATM 3835 FE   HEM     1      17.140   3.115  15.066  1.00 14.14          FE3+

const _Format_PDB_ATOM = FormatExpr(
    #          1         2         3         4         5         6         7         8
    # 12345678901234567890123456789012345678901234567890123456789012345678901234567890
    # >    <>   < >  <|> < |>  <|   >      <>      <>      <>    <>    <      >  <><><
    "{:<6}{:>5d} {:<4}{:>1}{:>3} {:>1}{:>4}{:>1}   {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}      {:<4}{:>2}{:>2}\n"
    )

# Models are numbered sequentially beginning with 1.
# Each MODEL must have a corresponding ENDMDL record.
# COLUMNS        DATA TYPE       CONTENTS
# --------------------------------------------------------------------------------
# 1 -  6       Record name    "MODEL "
# 11 - 14       Integer        Model serial number
# Example:
#          1         2         3         4         5         6         7         8
# 12345678901234567890123456789012345678901234567890123456789012345678901234567890
# MODEL        1
# ATOM      1  N   ALA     1      11.104   6.134  -6.504  1.00  0.00           N
# ATOM    294 2HG  GLU    18     -13.630  -3.769   0.160  1.00  0.00           H
# TER     295      GLU    18
# ENDMDL

const _Format_PDB_MODEL = FormatExpr(
    #          1
    # 12345678901234
    # MODEL        1
     "MODEL     {:>4}\n"
    )

# TER
# Indicates the end of a list of ATOM/HETATM records for a chain
# The TER records occur in the coordinate section of the entry, and indicate
# the last residue presented for each polypeptide and/or nucleic acid chain for
# which there are coordinates.
# The TER record has the same residue name, chain identifier, sequence number
# and insertion code as the terminal residue. The serial number of the TER
# record is one number greater than the serial number of the ATOM/HETATM
# preceding the TER.
# The residue name appearing on the TER record must be the same as the residue name
# of the immediately preceding ATOM or non-water HETATM record.
# COLUMNS         DATA TYPE         CONTENTS
# --------------------------------------------------------------------------------
#  1 -  6         Record name       "TER   "
#  7 - 11         Integer           Serial number
# 18 - 20         Residue name      Residue name
# 22              Character         Chain identifier
# 23 - 26         Integer           Residue sequence number
# 27              AChar             Insertion code
# Example:
#          1         2         3         4         5         6         7         8
# 12345678901234567890123456789012345678901234567890123456789012345678901234567890
# ATOM   4150  H   ALA A 431       8.674  16.036  12.858  1.00  0.00           H
# TER    4151      ALA A 431

# HETATM 1415  O2  BLE P   1      13.775  30.147  14.862  1.09 20.95           O
# TER    1416      BLE P   1

const _Format_PDB_TER = FormatExpr(
    # TER    4151      ALA A 431
    #          1         2         3         4         5         6         7         8
    # 12345678901234567890123456789012345678901234567890123456789012345678901234567890
    "TER   {:>5d}      {:>3} {:>1}{:>4}{:>1}\n"
    )

function print(io::IO, res::PDBResidue, format::Type{PDBFile}, atom_index::Int, serial_number::Int)
    number = match(r"(\d+)(\D?)", res.id.number)
    atomname = res.atoms[atom_index].atom
    printfmt(io, _Format_PDB_ATOM,
             res.id.group,
             serial_number,
             length(atomname) <= 3 ? string(" ", atomname) : atomname, # It works with NACCESS
             " ",
             res.id.name,
             res.id.chain,
             number[1],
             number[2],
             res.atoms[atom_index].coordinates.x,
             res.atoms[atom_index].coordinates.y,
             res.atoms[atom_index].coordinates.z,
             res.atoms[atom_index].occupancy,
             res.atoms[atom_index].B,
             " ",
             res.atoms[atom_index].element,
             " ")
    serial_number + 1
end

function print(io::IO, res::PDBResidue, format::Type{PDBFile}, start::Int=1)
    next = start
    for i in 1:length(res.atoms)
        next = print(io, res, format, i, next)
    end
    nothing
end

function print(io::IO, reslist::AbstractVector{PDBResidue}, format::Type{PDBFile}, start::Int=1)
    next = start

    use_model = length(unique(map(res->res.id.model, reslist))) > 1
    if use_model
        model = "START"
    end

    for resindex in 1:length(reslist)
        res = reslist[resindex]

        # MODEL
        if use_model
            if model != res.id.model
                if model != "START"
                   println(io, "ENDMDL")
                end
                printfmt(io, _Format_PDB_MODEL, res.id.model)
            end
            model = res.id.model
        end

        # TER

        # MIToS only prints TER for the ATOM group if the chain changes.
        # Some modified residues are annotated as HETATM in the middle of the ATOM chain:
        # TER can not be printed from ATOM to HETATM if the chain doesn’t change.
        if resindex > 1
            previous_res = reslist[ resindex - 1 ]
            if (previous_res.id.group == "ATOM") && (previous_res.id.chain != res.id.chain)
                number = match(r"(\d+)(\D?)", previous_res.id.number)
                printfmt(io, _Format_PDB_TER,
                    next,
                    previous_res.id.name,
                    previous_res.id.chain,
                    number[1],
                    number[2]
                    )
                next += 1
            end
        end

        # ATOM/HETATM
        for i in 1:length(res.atoms)
            next = print(io, res, format, i, next)
        end

    end

    if use_model
        println(io, "ENDMDL")
    end

    println(io, "END")
    nothing
end

print(reslist::AbstractVector{PDBResidue}, format::Type{PDBFile}) = print(STDOUT, reslist, format)
print(res::PDBResidue, format::Type{PDBFile}) = print(STDOUT, res, format)

@doc """
`print(io, res, format::Type{PDBFile})`
`print(res, format::Type{PDBFile})`

Print a `PDBResidue` or a vector of `PDBResidue`s in PDB format.
""" print
