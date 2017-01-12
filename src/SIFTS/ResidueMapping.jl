abstract DataBase

@auto_hash_equals immutable dbPDBe <: DataBase
    number::Int # Cross referenced residue number
    name::String # Cross referenced residue name
end

@inline _number_type(::Type{dbPDBe}) = Int

@auto_hash_equals immutable dbInterPro <: DataBase
    id::String
    number::String # Cross referenced residue number
    name::String # Cross referenced residue name
    evidence::String
end

@inline _number_type(::Type{dbInterPro}) = String

for ref_type in [:dbUniProt, :dbPfam, :dbNCBI]
    @eval begin

        @auto_hash_equals immutable $(ref_type) <: DataBase
            id::String # The cross reference database identifier
            number::Int # Cross referenced residue number
            name::String # Cross referenced residue name
        end

        @inline _number_type(::Type{$(ref_type)}) = Int
    end
end

for ref_type in [:dbPDB, :dbCATH, :dbSCOP]
    @eval begin

        @auto_hash_equals immutable $(ref_type) <: DataBase
            id::String
            number::String
            name::String
            chain::String
        end

        @inline _number_type(::Type{$(ref_type)}) = String
    end
end

"""Returns "" if the attributte is missing"""
function _get_attribute(elem::LightXML.XMLElement, attr::String)
    text = attribute(elem, attr)
    if text === nothing || text == "None"
        return("")
    else
        return(text)
    end
end

for ref_type in [:dbPDB, :dbCATH, :dbSCOP]
    @eval begin

        function (::Type{$(ref_type)})(map::LightXML.XMLElement)
            $(ref_type)(
                _get_attribute(map, "dbAccessionId"),
                _get_attribute(map, "dbResNum"),
                _get_attribute(map, "dbResName"),
                _get_attribute(map, "dbChainId")
            )
        end
    end
end

for ref_type in [:dbUniProt, :dbPfam, :dbNCBI]
    @eval begin

        function (::Type{$(ref_type)})(map::LightXML.XMLElement)
            $(ref_type)(
                _get_attribute(map, "dbAccessionId"),
                parse(Int, _get_attribute(map, "dbResNum")),
                _get_attribute(map, "dbResName")
            )
        end
    end
end

function (::Type{dbInterPro})(map::LightXML.XMLElement)
    dbInterPro(
        _get_attribute(map, "dbAccessionId"),
        _get_attribute(map, "dbResNum"),
        _get_attribute(map, "dbResName"),
        _get_attribute(map, "dbEvidence")
    )
end

function (::Type{dbPDBe})(map::LightXML.XMLElement)
      dbPDBe(
        parse(Int, _get_attribute(map, "dbResNum")),
        _get_attribute(map, "dbResName")
      )
end

@auto_hash_equals immutable SIFTSResidue
    PDBe::dbPDBe
    # crossRefDb
    UniProt::Nullable{dbUniProt}
    Pfam::Nullable{dbPfam}
    NCBI::Nullable{dbNCBI}
    InterPro::Array{dbInterPro, 1}
    PDB::Nullable{dbPDB}
    SCOP::Nullable{dbSCOP}
    CATH::Nullable{dbCATH}
    # residueDetail
    missing::Bool  # XML: <residueDetail dbSource="PDBe" property="Annotation" ...
    sscode::String # XML: <residueDetail dbSource="PDBe" property="codeSecondaryStructure"...
    ssname::String # XML: <residueDetail dbSource="PDBe" property="nameSecondaryStructure"...
end


# Getters
# -------

@inline _name(::Type{dbPDBe})    = "PDBe"
@inline _name(::Type{dbUniProt}) = "UniProt"
@inline _name(::Type{dbPfam})    = "Pfam"
@inline _name(::Type{dbNCBI})    = "NCBI"
@inline _name(::Type{dbInterPro})= "InterPro"
@inline _name(::Type{dbPDB})     = "PDB"
@inline _name(::Type{dbSCOP})    = "SCOP"
@inline _name(::Type{dbCATH})    = "CATH"

@inline Base.get(res::SIFTSResidue, db::Type{dbPDBe})    = res.PDBe
@inline Base.get(res::SIFTSResidue, db::Type{dbUniProt}) = res.UniProt
@inline Base.get(res::SIFTSResidue, db::Type{dbPfam})    = res.Pfam
@inline Base.get(res::SIFTSResidue, db::Type{dbNCBI})    = res.NCBI
@inline Base.get(res::SIFTSResidue, db::Type{dbInterPro})= res.InterPro
@inline Base.get(res::SIFTSResidue, db::Type{dbPDB})     = res.PDB
@inline Base.get(res::SIFTSResidue, db::Type{dbSCOP})    = res.SCOP
@inline Base.get(res::SIFTSResidue, db::Type{dbCATH})    = res.CATH

function Base.get{T<:Union{dbUniProt,dbPfam,dbNCBI,dbPDB,dbSCOP,dbCATH}}(res::SIFTSResidue,
                 db::Type{T}, field::Symbol, default::Union{String,Int})
    database = get(res, db)
    isnull(database) ? default : getfield(get(database), field)
end

function Base.get{T<:Union{dbUniProt,dbPfam,dbNCBI,dbPDB,dbSCOP,dbCATH}}(res::SIFTSResidue,
                 db::Type{T}, field::Symbol)
    database = get(res, db)
    S = fieldtype(T, field)
    isnull(database) ? Nullable{S}() : Nullable{S}(getfield(get(database),field))
end

# Print
# -----

function Base.show(io::IO, res::SIFTSResidue)
    if res.missing
        println(io, "SIFTSResidue (missing)")
    else
        println(io, "SIFTSResidue with secondary structure code (sscode): \"",
                res.sscode, "\" and name (ssname): \"", res.ssname, "\"")
    end
    println(io, "  PDBe:")
    println(io, "    number: ", res.PDBe.number)
    println(io, "    name: ", res.PDBe.name)
    for dbname in [:UniProt, :Pfam, :NCBI, :PDB, :SCOP, :CATH]
        dbfield = getfield(res, dbname)
        if !isnull(dbfield)
            println(io, "  ", dbname, " (Nullable) :")
            for f in fieldnames(get(dbfield))
                println(io, "    ", f, ": ",  getfield(get(dbfield), f))
            end
        end
    end
    println(io, "  InterPro: ",  res.InterPro)
end

# Creation
# --------

function (::Type{SIFTSResidue})(residue::LightXML.XMLElement, missing::Bool,
                                sscode::String, ssname::String)
    PDBe = dbPDBe(residue)
    UniProt = Nullable{dbUniProt}()
    Pfam = Nullable{dbPfam}()
    NCBI = Nullable{dbNCBI}()
    InterPro = dbInterPro[]
    PDB = Nullable{dbPDB}()
    SCOP = Nullable{dbSCOP}()
    CATH = Nullable{dbCATH}()
    for crossref in get_elements_by_tagname(residue, "crossRefDb")
        db = attribute(crossref, "dbSource")
        if db == "UniProt"
            UniProt = Nullable(dbUniProt(crossref))
        elseif db == "Pfam"
            Pfam = Nullable(dbPfam(crossref))
        elseif db == "NCBI"
            NCBI = Nullable(dbNCBI(crossref))
        elseif db == "InterPro"
            push!(InterPro, dbInterPro(crossref))
        elseif db == "PDB"
            PDB = Nullable(dbPDB(crossref))
        elseif db == "SCOP"
            SCOP = Nullable(dbSCOP(crossref))
        elseif db == "CATH"
            CATH = Nullable(dbCATH(crossref))
        else
            warn(string(db, " is not in the MIToS' DataBases."))
        end
    end
    SIFTSResidue(PDBe,
                 UniProt,
                 Pfam,
                 NCBI,
                 InterPro,
                 PDB,
                 SCOP,
                 CATH,
                 missing,
                 sscode,
                 ssname)
end

function (::Type{SIFTSResidue})(residue::LightXML.XMLElement)
    SIFTSResidue(residue, _get_details(residue)...)
end

# Mapping Functions
# =================

_is_All(::Any) = false
_is_All(::Type{All}) = true

_parse(::Type{Int}, str) = parse(Int, str)
_parse(::Type{String}, str) = ascii(str)
@inline _parse(::Type{String}, str::String) = str

"""
Parses a SIFTS XML file and returns a `Dict` between residue numbers of two `DataBase`s  with the given identifiers.
A `chain` could be specified (`All` by default).
If `missings` is `true` (default) all the residues are used, even if they haven’t coordinates in the PDB file.
"""
function siftsmapping{F, T}(filename::String,
                            db_from::Type{F},
                            id_from::String,
                            db_to::Type{T},
                            id_to::String;
                            chain::Union{Type{All},String} = All,
                            missings::Bool = true)
    mapping = Dict{_number_type(F), _number_type(T)}()
    xdoc = parse_file(filename)
    try
        for entity in _get_entities(xdoc)
            segments = _get_segments(entity)
            for segment in segments
                residues = _get_residues(segment)
                for residue in residues
                    in_chain = _is_All(chain)
                    key_data = _name(db_from) == "PDBe" ? Nullable(parse(Int,attribute(residue,"dbResNum"))) : Nullable{_number_type(F)}()
                    value_data = _name(db_to) == "PDBe" ? Nullable(parse(Int,attribute(residue,"dbResNum"))) : Nullable{_number_type(T)}()
                    if missings || !_is_missing(residue)
                        crossref = get_elements_by_tagname(residue, "crossRefDb")
                        for ref in crossref
                            source = attribute(ref, "dbSource")
                            if source == _name(db_from) && attribute(ref, "dbAccessionId") == id_from
                                key_data = Nullable(_parse(_number_type(F),attribute(ref,"dbResNum")))
                            end
                            if source == _name(db_to) && attribute(ref, "dbAccessionId") == id_to
                                value_data = Nullable(_parse(_number_type(T),attribute(ref,"dbResNum")))
                            end
                            if !in_chain && source == "PDB" # XML: <crossRefDb dbSource="PDB" ... dbChainId="E"/>
                                in_chain = attribute(ref, "dbChainId") == chain
                            end
                        end
                        if !isnull(key_data) && !isnull(value_data) && in_chain
                            key = get(key_data)
                            if haskey(mapping, key)
                                warn(string("$key is already in the mapping with the value ",
                                            mapping[key],". The value is replaced by ",
                                            get(value_data)))
                            end
                            mapping[key] = get(value_data)
                        end
                    end
                end
            end
        end
    finally
        free(xdoc)
    end
    sizehint!(mapping, length(mapping))
end

"""
`parse(document::LightXML.XMLDocument, ::Type{SIFTSXML}; chain=All, missings::Bool=true)`

Returns a `Vector{SIFTSResidue}` parsed from a `SIFTSXML` file.
By default, parses all the `chain`s and includes `missings` residues.
"""
function Base.parse(document::LightXML.XMLDocument, ::Type{SIFTSXML};
                    chain::Union{Type{All},String}=All, missings::Bool = true)
    vector = SIFTSResidue[]
    for entity in _get_entities(document)
        for segment in _get_segments(entity)
            residues = _get_residues(segment)
            for residue in residues
                missing, sscode, ssname = _get_details(residue)
                if missings || !missing
                    sifts_res = SIFTSResidue(residue, missing, sscode, ssname)
                    if _is_All(chain) || (!isnull(sifts_res.PDB) && get(sifts_res.PDB).chain == chain)
                        push!(vector, sifts_res)
                    end
                end
            end
        end
    end
    vector
end

# Find SIFTSResidue
# -----------------

for F in (:find, :filter!, :filter)
    @eval begin
        function Base.$(F){T<:DataBase}(f::Function,list::AbstractVector{SIFTSResidue},db::Type{T})
            $(F)(list) do res
                database = get(res, db)
                if !isnull(database)
                    f(get(database))
                end
            end
        end
    end
end
