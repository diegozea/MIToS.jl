immutable PDBFile <: Format end

"""
Reads a text file of a PDB entry.
Returns a list of PDBResidue (view MIToS.PDB.PDBResidues).
Setting `chain`, `model`, `group`, `atomname` and `onlyheavy` values
can be used to select of a subset of all residues. Group can be ATOM
or HETATM. If not set, all residues are returned.
If the keyword argument `occupancyfilter` (default: `false`) is `true`,
only the atoms with the best occupancy are returned.
"""
function parse(io::Union{IO, ASCIIString}, ::Type{PDBFile}; chain::ASCIIString = "all",
               model::ASCIIString = "all", group::ASCIIString = "all", atomname::ASCIIString="all",
               onlyheavy::Bool=false, occupancyfilter::Bool=false)
    residue_dict = OrderedDict{PDBResidueIdentifier, Vector{PDBAtom}}()
    atom_model = 0
    for line in eachline(io)
        line_id = replace(line[1:6], ' ', "")
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
    for res in reslist
        for i in 1:length(res.atoms)
            next = print(io, res, format, i, next)
        end
    end
    nothing
end

print(reslist::AbstractVector{PDBResidue}, format::Type{PDBFile}) = print(STDOUT, reslist, format)
print(res::PDBResidue, format::Type{PDBFile}) = print(STDOUT, res, format)
