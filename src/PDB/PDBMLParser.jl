"""
`PDBML <: FileFormat`

Protein Data Bank Markup Language (PDBML), a representation of PDB data in XML format.
"""
struct PDBML <: FileFormat end

function _get_text(elem, name)
    sub = find_element(elem, name)
    if sub !== nothing
        return(content(sub))
    end
    throw(ErrorException("There is not $name for $elem"))
end

function _get_ins_code(elem)
    sub = find_element(elem, "pdbx_PDB_ins_code")
    if sub !== nothing
        return(content(sub))
    else
        return("")
    end
end

function _get_atom_iterator(document::LightXML.XMLDocument)
    pdbroot = root(document)
    child_elements(get_elements_by_tagname(pdbroot, "atom_siteCategory")[1])
end

"Used for parsing a PDB file into `Vector{PDBResidue}`"
function _generate_residues(residue_dict::OrderedDict{PDBResidueIdentifier, Vector{PDBAtom}},
                            occupancyfilter::Bool=false)
    if occupancyfilter
        return( PDBResidue[ PDBResidue(k, bestoccupancy(v)) for (k,v) in residue_dict ] )
    else
        return( PDBResidue[ PDBResidue(k, v) for (k,v) in residue_dict ] )
    end
end

"""
`parse(pdbml, ::Type{PDBML}; chain=All, model=All, group=All, atomname=All, onlyheavy=false, label=true, occupancyfilter=false)`

Reads a `LightXML.XMLDocument` representing a pdb file.
Returns a list of `PDBResidue`s (view `MIToS.PDB.PDBResidues`).
Setting `chain`, `model`, `group`, `atomname` and `onlyheavy` values can be used to select
of a subset of all residues. If not set, all residues are returned. If the keyword
argument `label` (default: `true`) is `false`,the **auth_** attributes will be use instead
of the **label_** attributes for `chain`, `atom` and residue `name` fields. The **auth_**
attributes are alternatives provided by an author in order to match the
identification/values used in the publication that describes the structure. If the
keyword argument `occupancyfilter` (default: `false`) is `true`, only the atoms with the
best occupancy are returned.
"""
function Base.parse(pdbml::LightXML.XMLDocument, ::Type{PDBML};
                    chain::Union{String,Type{All}} = All,
                    model::Union{String,Type{All}} = All,
                    group::Union{String,Type{All}} = All,
                    atomname::Union{String,Type{All}} = All,
                    onlyheavy::Bool=false,
                    label::Bool=true,
                    occupancyfilter::Bool=false)

    residues = Vector{PDBResidue}()

    prefix = label ? "label" : "auth"
    chain_attribute = string(prefix, "_asym_id")
    atom_attribute = string(prefix, "_atom_id")
    comp_attribute = string(prefix, "_comp_id")

    residue_id = PDBResidueIdentifier("", "", "", "", "", "")
    atoms = _get_atom_iterator(pdbml)
    for atom in atoms

        atom_name = _get_text(atom, atom_attribute)
        if !_is(atom_name, atomname)
            continue
        end

        element = _get_text(atom, "type_symbol")
        if onlyheavy && element=="H"
            continue
        end

        atom_group = _get_text(atom, "group_PDB")
        if !_is(atom_group, group)
            continue
        end

        atom_chain = _get_text(atom, chain_attribute)
        if !_is(atom_chain, chain)
            continue
        end

        atom_model = _get_text(atom, "pdbx_PDB_model_num")
        if !_is(atom_model, model)
            continue
        end

        PDBe_number = _get_text(atom, "label_seq_id")

        #  Residue_No  _atom_site.auth_seq_id
        #  Ins_Code    _atom_site.pdbx_PDB_ins_code
        PDB_number = string(_get_text(atom, "auth_seq_id"), _get_ins_code(atom))
        name = _get_text(atom, comp_attribute)

        if (residue_id.PDBe_number != PDBe_number) ||
                (residue_id.number != PDB_number) ||
                (residue_id.name != name) ||
                (residue_id.chain != atom_chain) ||
                (residue_id.group != atom_group) ||
                (residue_id.model != atom_model)

            n_res = length(residues)
            if occupancyfilter && n_res > 0
                residues[n_res].atoms = bestoccupancy(residues[n_res].atoms)
            end

            residue_id = PDBResidueIdentifier(PDBe_number,
                                              PDB_number,
                                              name,
                                              atom_group,
                                              atom_model,
                                              atom_chain)
            push!(residues, PDBResidue(residue_id, Vector{PDBAtom}()))
        end

        x = parse(Float64, _get_text(atom, "Cartn_x"))
        y = parse(Float64, _get_text(atom, "Cartn_y"))
        z = parse(Float64, _get_text(atom, "Cartn_z"))
        occupancy = parse(Float64, _get_text(atom, "occupancy"))
        B = _get_text(atom, "B_iso_or_equiv")

        push!(residues[end].atoms, PDBAtom(Coordinates(x,y,z),
                                           atom_name,
                                           element,
                                           occupancy,
                                           B))
    end

    if occupancyfilter
        residues[end].atoms = bestoccupancy(residues[end].atoms)
    end

    residues
end

# function Base.parse(pdbml::LightXML.XMLDocument, ::Type{PDBML};
#                     chain::Union{String,Type{All}} = All,
#                     model::Union{String,Type{All}} = All,
#                     group::Union{String,Type{All}} = All,
#                     atomname::Union{String,Type{All}} = All,
#                     onlyheavy::Bool=false,
#                     label::Bool=true,
#                     occupancyfilter::Bool=false)
#
#     residue_dict = OrderedDict{PDBResidueIdentifier, Vector{PDBAtom}}()
#
#     prefix = label ? "label" : "auth"
#     chain_attribute = string(prefix, "_asym_id")
#     atom_attribute = string(prefix, "_atom_id")
#     comp_attribute = string(prefix, "_comp_id")
#
#     atoms = _get_atom_iterator(pdbml)
#     for atom in atoms
#
#         atom_group = _get_text(atom, "group_PDB")
#         atom_model = _get_text(atom, "pdbx_PDB_model_num")
#         atom_chain = _get_text(atom, chain_attribute)
#         atom_name = _get_text(atom, atom_attribute)
#         element = _get_text(atom, "type_symbol")
#
#         if  _is(atom_group, group) && _is(atom_chain,chain) &&
#                 _is(atom_model, model) && _is(atom_name,atomname) &&
#                 (!onlyheavy || element!="H")
#
#             PDBe_number = _get_text(atom, "label_seq_id")
#
#             #  Residue_No  _atom_site.auth_seq_id
#             #  Ins_Code    _atom_site.pdbx_PDB_ins_code
#             PDB_number = string(_get_text(atom, "auth_seq_id"), _get_ins_code(atom))
#             name = _get_text(atom, comp_attribute)
#             x = parse(Float64, _get_text(atom, "Cartn_x"))
#             y = parse(Float64, _get_text(atom, "Cartn_y"))
#             z = parse(Float64, _get_text(atom, "Cartn_z"))
#             occupancy = parse(Float64, _get_text(atom, "occupancy"))
#             B = _get_text(atom, "B_iso_or_equiv")
#
#             residue_id = PDBResidueIdentifier(PDBe_number, PDB_number, name, atom_group, atom_model, atom_chain)
#             atom_data  = PDBAtom(Coordinates(x,y,z), atom_name, element, occupancy, B)
#
#             value = get!(residue_dict, residue_id, PDBAtom[])
#             push!(value, atom_data)
#         end
#     end
#     _generate_residues(residue_dict, occupancyfilter)
# end

# Download PDB
# ============

function _inputnameforgzip(outfile)
    if endswith(outfile, ".gz")
        return(outfile)
    end
    string(outfile, ".gz")
end

_file_extension(format::Type{PDBML}) = ".xml.gz"
_file_extension(format::Type{PDBFile}) = ".pdb.gz"

"""
It downloads a gzipped PDB file from PDB database.
It requires a four character `pdbcode`.
Its default `format` is `PDBML` (PDB XML) and It uses the
`baseurl` "http://www.rcsb.org/pdb/files/".
`filename` is the path/name of the output file.
This function calls `MIToS.Utils.download_file` that calls `HTTP.open`.
You can use keyword arguments from `HTTP.request` (e.g. `redirect`).
Use the `headers` keyword argument to pass a `Dict{String, String}`
with the header information.
"""
function downloadpdb(pdbcode::String;
                     format::Type{T} = PDBML,
                     filename::String= uppercase(pdbcode)*_file_extension(format),
                     baseurl::String = "http://www.rcsb.org/pdb/files/",
                     headers::Dict{String, String} = Dict{String, String}(),
                     kargs...) where T <: FileFormat
    if check_pdbcode(pdbcode)
        pdbfilename = uppercase(pdbcode) * _file_extension(format)
        filename = _inputnameforgzip(filename)
        sepchar = endswith(baseurl,"/") ? "" : "/";
        download_file(string(baseurl,sepchar,pdbfilename), filename; headers=headers, kargs...)
    else
        throw(ErrorException("$pdbcode is not a correct PDB code"))
    end
    filename
end


# RESTful PDB interface
# =====================

function downloadpdbheader(pdbcode::String; filename::String=tempname()*".xml", kargs...)
    pdbcode = lowercase(pdbcode)
    @assert endswith(filename,".xml") "filename must end with the xml extension: .xml"
    if check_pdbcode(pdbcode)
        query = "http://www.rcsb.org/pdb/rest/describePDB?structureId=" * pdbcode
        download_file(query, filename; kargs...)
    else
        throw(ErrorException("$pdbcode is not a correct PDB code"))
    end
    filename
end

struct PDBMLHeader <: FileFormat end

function Base.parse(document::LightXML.XMLDocument, ::Type{PDBMLHeader})
    xmlroot = root(document)
    description = get_elements_by_tagname(xmlroot, "PDB")[1]
    attributes_dict(description)
end

"""
Access general information about a PDB entry (e.g., Header information) using the
RESTful interface of the PDB database (describePDB). Returns a Dict for the four
character `pdbcode`.
"""
function getpdbdescription(pdbcode::String; kargs...)
    filename = downloadpdbheader(pdbcode)
    try
        read(filename, PDBMLHeader)
    finally
        rm(filename)
    end
end
