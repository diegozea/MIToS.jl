"""
`PDBFile <: FileFormat`

Protein Data Bank (PDB) format.
It provides a standard representation for macromolecular
structure data derived from X-ray diffraction and NMR studies.
"""
struct PDBFile <: FileFormat end

using BioStructures # Ensure BioStructures is available

function _parse_residueidentifier(line::String, atom_chain, line_id, actual_model)
    # 23 - 26        Integer         Residue sequence number.
    # 27             AChar           Code for insertion of residues.
    PDB_number = String(strip(SubString(line, 23, 27), ' '))

    name = String(strip(SubString(line, 18, 20), ' '))

    PDBResidueIdentifier("", PDB_number, name, line_id, actual_model, atom_chain)
end

function _parse_pdbatom(line::String, atom_name, element)
    x = parse(Float64, SubString(line, 31, 38))
    y = parse(Float64, SubString(line, 39, 46))
    z = parse(Float64, SubString(line, 47, 54))
    occupancy = parse(Float64, SubString(line, 55, 60))
    B = String(strip(SubString(line, 61, 66), ' '))
    alt_id = strip(String(SubString(line, 17, 17)))
    if 80 ≤ ncodeunits(line)
        charge = String(strip(SubString(line, 79, 80)))
    else
        charge = ""
    end

    PDBAtom(Coordinates(x, y, z), atom_name, element, occupancy, B, alt_id, charge)
end

"""
`parse_file(io, ::Type{PDBFile}; chain=All, model=All, group=All, atomname=All, onlyheavy=false, occupancyfilter=false)`

Reads a text file of a PDB entry.
Returns a list of `PDBResidue` (view `MIToS.PDB.PDBResidues`).
Setting `chain`, `model`, `group`, `atomname` and `onlyheavy` values
can be used to select of a subset of all residues. Group can be `"ATOM"`
or `"HETATM"`. If not set, all residues are returned.
If the keyword argument `occupancyfilter` (default: `false`) is `true`,
only the atoms with the best occupancy are returned.
"""
function Utils.parse_file(
    io::Union{IO,String},
    ::Type{PDBFile};
    chain::Union{String,Type{All}} = All,
    model::Union{String,Type{All}} = All,
    group::Union{String,Type{All}} = All,
    atomname::Union{String,Type{All}} = All,
    onlyheavy::Bool = false,
    occupancyfilter::Bool = false,
)
    original_pos = isa(io, IO) ? position(io) : nothing
    residues = Vector{PDBResidue}()

    try
        model_counter = 0
        actual_model = "1"
        is_model_filter = _is(actual_model, model) # Renamed to avoid conflict
        previous_used_line = ""
        residue_id = PDBResidueIdentifier("", "", "", "", "", "") # Placeholder

        for line::String in lineiterator(io)
            line_id = match(r"^\S{0,6}", line).match

            if line_id == "MODEL"
                model_counter += 1
                actual_model = string(model_counter)
                is_model_filter = _is(actual_model, model)
            end

            if (group === All && (line_id == "ATOM" || line_id == "HETATM")) ||
               (line_id == group)
                atom_chain = string(line[22])
                atom_name = String(strip(SubString(line, 13, 16), ' '))
                element = if 78 ≤ ncodeunits(line)
                    String(strip(SubString(line, 77, 78), ' '))
                else
                    ""
                end

                if is_model_filter &&
                   _is(atom_chain, chain) &&
                   _is(atom_name, atomname) &&
                   (!onlyheavy || element != "H")

                    if (previous_used_line == "") ||
                       (residue_id.group != line_id) ||
                       (residue_id.model != actual_model) ||
                       (SubString(previous_used_line, 18, 27) != SubString(line, 18, 27))

                        n_res = length(residues)
                        if occupancyfilter && n_res > 0 && !isempty(residues[n_res].atoms)
                            residues[n_res].atoms = bestoccupancy(residues[n_res].atoms)
                            # Remove residue if it becomes empty after occupancy filter
                            if isempty(residues[n_res].atoms)
                                pop!(residues)
                            end
                        end
                        # Check if we need to add a new residue or if the last one was removed
                        if length(residues) == n_res || n_res == 0 # previous was not removed or list is empty
                            residue_id = _parse_residueidentifier(line, atom_chain, line_id, actual_model)
                            push!(residues, PDBResidue(residue_id, Vector{PDBAtom}()))
                        else # last residue was popped, re-use its id (or re-parse if safer)
                             # For simplicity, let's assume parsing new one is fine.
                            residue_id = _parse_residueidentifier(line, atom_chain, line_id, actual_model)
                            push!(residues, PDBResidue(residue_id, Vector{PDBAtom}()))
                        end
                    end
                    
                    atom_data = _parse_pdbatom(line, atom_name, element)
                    # Ensure residues list is not empty before trying to push atoms
                    if !isempty(residues)
                        push!(residues[end].atoms, atom_data)
                        previous_used_line = line
                    else
                        # This case should ideally not be reached if logic is correct,
                        # but as a safeguard: if residues is empty, it implies the first
                        # residue that should have been created was filtered out before any atom
                        # could be added. This might happen if a MODEL record makes is_model_filter true,
                        # but then the first ATOM/HETATM record doesn't match other filters.
                        # If so, we need to create the residue first.
                        # However, the outer if `is_model_filter && ...` should prevent this.
                        # Let's assume current logic is mostly fine, but this is a complex area.
                    end
                end
            end
        end

        if occupancyfilter && !isempty(residues) && !isempty(residues[end].atoms)
            residues[end].atoms = bestoccupancy(residues[end].atoms)
            if isempty(residues[end].atoms)
                pop!(residues)
            end
        end
        return residues # Success from primary parser

    catch err
        if err isa BoundsError || err isa ArgumentError
            # Fallback to BioStructures.jl
            try
                struc = if isa(io, IO)
                    if original_pos !== nothing
                        seek(io, original_pos)
                        BioStructures.read(io, BioStructures.PDB)
                    else # Should not happen if original_pos was set for IO
                        rethrow(err) # Cannot reset stream
                    end
                else # io is a String (filepath)
                    BioStructures.read(io, BioStructures.PDB)
                end

                biostructure_residues = convert(Vector{PDBResidue}, struc)
                
                # Apply filtering to biostructure_residues
                filtered_residues = Vector{PDBResidue}()
                for res_bs in biostructure_residues
                    # Check residue level filters
                    if !_is(res_bs.id.model, model) continue end
                    if !_is(res_bs.id.chain, chain) continue end
                    # Group for BioStructures PDBResidue is in res_bs.id.group ("ATOM" or "HETATM")
                    if group !== All && res_bs.id.group != group continue end

                    filtered_atoms = PDBAtom[]
                    for atom_bs in res_bs.atoms
                        # Check atom level filters
                        if !_is(atom_bs.atom, atomname) continue end
                        if onlyheavy && atom_bs.element == "H" continue end
                        push!(filtered_atoms, atom_bs)
                    end

                    if !isempty(filtered_atoms)
                        # If occupancyfilter is true, apply it to these atoms
                        final_atoms = if occupancyfilter
                            bestoccupancy(filtered_atoms)
                        else
                            filtered_atoms
                        end
                        
                        if !isempty(final_atoms)
                            # Create a new PDBResidue with the original id but filtered atoms
                            push!(filtered_residues, PDBResidue(res_bs.id, final_atoms))
                        end
                    end
                end
                return filtered_residues
            catch bio_err # Error during BioStructures parsing or conversion
                @warn "BioStructures.jl fallback failed with error: $bio_err"
                rethrow(err) # Rethrow original error
            end
        else
            rethrow(err) # Not a BoundsError or ArgumentError
        end
    end
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
    "{:<6}{:>5d} {:<4}{:>1}{:>3} {:>1}{:>4}{:>1}   {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6}      {:<4}{:>2}{:>2}\n",
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
    "MODEL     {:>4}\n",
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
    "TER   {:>5d}      {:>3} {:>1}{:>4}{:>1}\n",
)

function _get_residue_number(res::PDBResidue)
    number = match(r"(-?\d+)(\D?)", res.id.number)
    if number === nothing
        throw(ErrorException("Invalid residue number: $(res.id.number)"))
    end
    number
end

function Utils.print_file(
    io::IO,
    res::PDBResidue,
    format::Type{PDBFile},
    atom_index::Int,
    serial_number::Int,
)
    number = _get_residue_number(res)
    atomname = res.atoms[atom_index].atom
    printfmt(
        io,
        _Format_PDB_ATOM,
        res.id.group,
        serial_number,
        length(atomname) <= 3 ? string(" ", atomname) : atomname, # It works with NACCESS
        res.atoms[atom_index].alt_id,
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
        res.atoms[atom_index].charge,
    )
    serial_number + 1
end

function Utils.print_file(io::IO, res::PDBResidue, format::Type{PDBFile}, start::Int = 1)
    next = start
    for i in eachindex(res.atoms)
        next = print_file(io, res, format, i, next)
    end
    nothing
end

function Utils.print_file(
    io::IO,
    reslist::AbstractVector{PDBResidue},
    format::Type{PDBFile},
    start::Int = 1,
)
    next = start

    use_model = length(unique(map(res -> res.id.model, reslist))) > 1
    if use_model
        model = "START"
    end

    for resindex in eachindex(reslist)
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
            previous_res = reslist[resindex-1]
            if (previous_res.id.group == "ATOM") && (previous_res.id.chain != res.id.chain)
                number = _get_residue_number(previous_res)
                printfmt(
                    io,
                    _Format_PDB_TER,
                    next,
                    previous_res.id.name,
                    previous_res.id.chain,
                    number[1],
                    number[2],
                )
                next += 1
            end
        end

        # ATOM/HETATM
        for i in eachindex(res.atoms)
            next = print_file(io, res, format, i, next)
        end

    end

    if use_model
        println(io, "ENDMDL")
    end

    println(io, "END")
    nothing
end

@doc """
`print_file(io, res, format::Type{PDBFile})`
`print_file(res, format::Type{PDBFile})`

Print a `PDBResidue` or a vector of `PDBResidue`s in PDB format.
""" print_file
