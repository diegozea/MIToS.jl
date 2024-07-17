abstract type DataBase end

"""
`dbPDBe` stores the residue `number` and `name` in PDBe as strings.
"""
@auto_hash_equals struct dbPDBe <: DataBase
    number::String # Cross referenced residue number
    name::String # Cross referenced residue name
end

"""
`dbInterPro` stores the residue `id`, `number`, `name` and `evidence` in InterPro as strings.
"""
@auto_hash_equals struct dbInterPro <: DataBase
    id::String
    number::String # Cross referenced residue number
    name::String # Cross referenced residue name
    evidence::String
end

"""
`dbEnsembl` stores the residue (gene) accession `id`, the `transcript`,
`translation` and `exon` ids in Ensembl as strings, together with the residue
`number` and `name` using the UniProt coordinates.
"""
@auto_hash_equals struct dbEnsembl <: DataBase
    id::String # (gene) accession id
    number::String # Cross referenced residue number
    name::String # Cross referenced residue name
    transcript::String
    translation::String
    exon::String
end

for ref_type in [:dbUniProt, :dbPfam, :dbNCBI]
    @eval begin

        @auto_hash_equals struct $(ref_type) <: DataBase
            id::String # The cross reference database identifier
            number::String # Cross referenced residue number
            name::String # Cross referenced residue name
        end

    end
end

@doc """
`dbUniProt` stores the residue `id`, `number` and `name` in UniProt as strings.
""" dbUniProt

@doc """
`dbPfam` stores the residue `id`, `number` and `name` in Pfam as strings.
""" dbPfam

@doc """
`dbNCBI` stores the residue `id`, `number` and `name` in NCBI as strings.
""" dbNCBI

for ref_type in [:dbPDB, :dbCATH, :dbSCOP, :dbSCOP2, :dbSCOP2B]
    @eval begin

        @auto_hash_equals struct $(ref_type) <: DataBase
            id::String
            number::String
            name::String
            chain::String
        end

    end
end

@doc """
`dbPDB` stores the residue `id`, `number`, `name` and `chain` in PDB as strings.
""" dbPDB

@doc """
`dbCATH` stores the residue `id`, `number`, `name` and `chain` in CATH as 
strings.
""" dbCATH

@doc """
`dbSCOP` stores the residue `id`, `number`, `name` and `chain` in SCOP as 
strings.
""" dbSCOP

@doc """
`dbSCOP2` stores the residue `id`, `number`, `name` and `chain` in SCOP2 as 
strings.
""" dbSCOP2

@doc """
`dbSCOP2B` stores the residue `id`, `number`, `name` and `chain` in SCOP2B as 
strings. *SCOP2B* is expansion of *SCOP2* domain annotations at superfamily 
level to every *PDB* with same *UniProt* accession having at least 80% *SCOP2* 
domain coverage.
""" dbSCOP2B

"""
Returns "" if the attributte is missing
"""
function _get_attribute(elem::LightXML.XMLElement, attr::String)
    text = LightXML.attribute(elem, attr)
    if text === nothing || text == "None"
        return ("")
    else
        return (text)
    end
end

"""
Returns `missing` if the attributte is missing
"""
function _get_nullable_attribute(
    elem::LightXML.XMLElement,
    attr::String,
)::Union{String,Missing}
    text = LightXML.attribute(elem, attr)
    (text === nothing || text == "None") ? missing : text
end

for ref_type in [:dbPDB, :dbCATH, :dbSCOP, :dbSCOP2, :dbSCOP2B]
    @eval begin

        function $(ref_type)(map::LightXML.XMLElement)
            $(ref_type)(
                _get_attribute(map, "dbAccessionId"),
                _get_attribute(map, "dbResNum"),
                _get_attribute(map, "dbResName"),
                _get_attribute(map, "dbChainId"),
            )
        end
    end
end

for ref_type in [:dbUniProt, :dbPfam, :dbNCBI]
    @eval begin

        function $(ref_type)(map::LightXML.XMLElement)
            $(ref_type)(
                _get_attribute(map, "dbAccessionId"),
                _get_attribute(map, "dbResNum"),
                _get_attribute(map, "dbResName"),
            )
        end
    end
end


function dbEnsembl(map::LightXML.XMLElement)
    dbEnsembl(
        _get_attribute(map, "dbAccessionId"),
        _get_attribute(map, "dbResNum"),
        _get_attribute(map, "dbResName"),
        _get_attribute(map, "dbTranscriptId"),
        _get_attribute(map, "dbTranslationId"),
        _get_attribute(map, "dbExonId"),
    )
end


function dbInterPro(map::LightXML.XMLElement)
    dbInterPro(
        _get_attribute(map, "dbAccessionId"),
        _get_attribute(map, "dbResNum"),
        _get_attribute(map, "dbResName"),
        _get_attribute(map, "dbEvidence"),
    )
end

function dbPDBe(map::LightXML.XMLElement)
    dbPDBe(_get_attribute(map, "dbResNum"), _get_attribute(map, "dbResName"))
end

"""
A `SIFTSResidue` object stores the SIFTS residue level mapping for a residue. It has the
following fields that you can access at any moment for query purposes:

    - `PDBe` : A `dbPDBe` object, it's present in all the `SIFTSResidue`s.
    - `UniProt` : A `dbUniProt` object or `missing`.
    - `Pfam` : A `dbPfam` object or `missing`.
    - `NCBI` : A `dbNCBI` object or `missing`.
    - `InterPro` : An array of `dbInterPro` objects.
    - `PDB` : A `dbPDB` object or `missing`.
    - `SCOP` : A `dbSCOP` object or `missing`.
    - `SCOP2` : An array of `dbSCOP2` objects.
    - `SCOP2B` : A `dbSCOP2B` object or `missing`.
    - `CATH` : A `dbCATH` object or `missing`.
    - `Ensembl` : An array of `dbEnsembl` objects.
    - `missing` : It's `true` if the residue is missing, i.e. not observed, in the structure.
    - `sscode` : A string with the secondary structure code of the residue.
    - `ssname` : A string with the secondary structure name of the residue.
"""
@auto_hash_equals struct SIFTSResidue
    PDBe::dbPDBe
    # crossRefDb
    UniProt::Union{dbUniProt,Missing}
    Pfam::Union{dbPfam,Missing}
    NCBI::Union{dbNCBI,Missing}
    InterPro::Array{dbInterPro,1}
    PDB::Union{dbPDB,Missing}
    SCOP::Union{dbSCOP,Missing}
    SCOP2::Array{dbSCOP2,1}
    SCOP2B::Union{dbSCOP2B,Missing}
    CATH::Union{dbCATH,Missing}
    Ensembl::Array{dbEnsembl,1}
    # residueDetail
    missing::Bool  # XML: <residueDetail dbSource="PDBe" property="Annotation" ...
    sscode::String # XML: <residueDetail dbSource="PDBe" property="codeSecondaryStructure"...
    ssname::String # XML: <residueDetail dbSource="PDBe" property="nameSecondaryStructure"...
end


# Getters
# -------

@inline _name(::Type{dbPDBe}) = "PDBe"
@inline _name(::Type{dbUniProt}) = "UniProt"
@inline _name(::Type{dbPfam}) = "Pfam"
@inline _name(::Type{dbNCBI}) = "NCBI"
@inline _name(::Type{dbInterPro}) = "InterPro"
@inline _name(::Type{dbPDB}) = "PDB"
@inline _name(::Type{dbSCOP}) = "SCOP"
@inline _name(::Type{dbSCOP2}) = "SCOP2"
@inline _name(::Type{dbSCOP2B}) = "SCOP2B"
@inline _name(::Type{dbCATH}) = "CATH"
@inline _name(::Type{dbEnsembl}) = "Ensembl"

@inline Base.get(res::SIFTSResidue, db::Type{dbPDBe}) = res.PDBe
@inline Base.get(res::SIFTSResidue, db::Type{dbUniProt}) = res.UniProt
@inline Base.get(res::SIFTSResidue, db::Type{dbPfam}) = res.Pfam
@inline Base.get(res::SIFTSResidue, db::Type{dbNCBI}) = res.NCBI
@inline Base.get(res::SIFTSResidue, db::Type{dbInterPro}) = res.InterPro
@inline Base.get(res::SIFTSResidue, db::Type{dbPDB}) = res.PDB
@inline Base.get(res::SIFTSResidue, db::Type{dbSCOP}) = res.SCOP
@inline Base.get(res::SIFTSResidue, db::Type{dbSCOP2}) = res.SCOP2
@inline Base.get(res::SIFTSResidue, db::Type{dbSCOP2B}) = res.SCOP2B
@inline Base.get(res::SIFTSResidue, db::Type{dbCATH}) = res.CATH
@inline Base.get(res::SIFTSResidue, db::Type{dbEnsembl}) = res.Ensembl

function Base.get(
    res::SIFTSResidue,
    db::Type{T},
    field::Symbol,
    default::Union{String,Missing} = missing,
) where {T<:Union{dbUniProt,dbPfam,dbNCBI,dbPDB,dbSCOP,dbSCOP2B,dbCATH}}
    database = get(res, db)
    ismissing(database) ? default : getfield(database, field)
end

# Print
# -----

function Base.show(io::IO, res::SIFTSResidue)
    if res.missing
        println(io, "SIFTSResidue (missing)")
    else
        println(
            io,
            "SIFTSResidue with secondary structure code (sscode): \"",
            res.sscode,
            "\" and name (ssname): \"",
            res.ssname,
            "\"",
        )
    end
    println(io, "  PDBe:")
    println(io, "    number: ", res.PDBe.number)
    println(io, "    name: ", res.PDBe.name)
    for dbname in [:UniProt, :Pfam, :NCBI, :PDB, :SCOP, :SCOP2B, :CATH]
        dbfield = getfield(res, dbname)
        if !ismissing(dbfield)
            println(io, "  ", dbname, " :")
            for f in fieldnames(typeof(dbfield))
                println(io, "    ", f, ": ", getfield(dbfield, f))
            end
        end
    end
    length(res.SCOP2) > 0 && println(io, "  SCOP2: ", res.SCOP2)
    length(res.InterPro) > 0 && println(io, "  InterPro: ", res.InterPro)
    length(res.Ensembl) > 0 && println(io, "  Ensembl: ", res.Ensembl)
end

# Creation
# --------

function SIFTSResidue(
    residue::LightXML.XMLElement,
    missing_residue::Bool,
    sscode::String,
    ssname::String,
)
    PDBe = dbPDBe(residue)
    UniProt = missing
    Pfam = missing
    NCBI = missing
    InterPro = dbInterPro[]
    PDB = missing
    SCOP = missing
    SCOP2 = dbSCOP2[]
    SCOP2B = missing
    CATH = missing
    Ensembl = dbEnsembl[]
    for crossref in LightXML.get_elements_by_tagname(residue, "crossRefDb")
        db = LightXML.attribute(crossref, "dbSource")
        if db == "UniProt"
            UniProt = dbUniProt(crossref)
        elseif db == "Pfam"
            Pfam = dbPfam(crossref)
        elseif db == "NCBI"
            NCBI = dbNCBI(crossref)
        elseif db == "InterPro"
            push!(InterPro, dbInterPro(crossref))
        elseif db == "PDB"
            PDB = dbPDB(crossref)
        elseif db == "SCOP"
            SCOP = dbSCOP(crossref)
        elseif db == "SCOP2"
            push!(SCOP2, dbSCOP2(crossref))
        elseif db == "SCOP2B"
            SCOP2B = dbSCOP2B(crossref)
        elseif db == "CATH"
            CATH = dbCATH(crossref)
        elseif db == "Ensembl"
            push!(Ensembl, dbEnsembl(crossref))
        else
            @warn(string(db, " is not in the MIToS' DataBases."))
        end
    end
    SIFTSResidue(
        PDBe,
        UniProt,
        Pfam,
        NCBI,
        InterPro,
        PDB,
        SCOP,
        SCOP2,
        SCOP2B,
        CATH,
        Ensembl,
        missing_residue,
        sscode,
        ssname,
    )
end

function SIFTSResidue(residue::LightXML.XMLElement)
    SIFTSResidue(residue, _get_details(residue)...)
end

# Mapping Functions
# =================

_is_All(::Any) = false
_is_All(::Type{All}) = true

"""
Parses a SIFTS XML file and returns a `OrderedDict` between residue numbers of
two `DataBase`s with the given identifiers. A `chain` could be specified
(`All` by default). If `missings` is `true` (default) all the residues are
used, even if they haven’t coordinates in the PDB file.
"""
function siftsmapping(
    filename::String,
    db_from::Type{F},
    id_from::String,
    db_to::Type{T},
    id_to::String;
    chain::Union{Type{All},String} = All,
    missings::Bool = true,
) where {F,T}
    mapping = OrderedDict{String,String}()
    xdoc = LightXML.parse_file(filename)
    try
        for entity in _get_entities(xdoc)
            segments = _get_segments(entity)
            for segment in segments
                residues = _get_residues(segment)
                for residue in residues
                    in_chain = _is_All(chain)
                    key_data =
                        _name(db_from) == "PDBe" ? LightXML.attribute(residue, "dbResNum") :
                        missing
                    value_data =
                        _name(db_to) == "PDBe" ? LightXML.attribute(residue, "dbResNum") :
                        missing
                    if missings || !_is_missing(residue)
                        crossref = LightXML.get_elements_by_tagname(residue, "crossRefDb")
                        for ref in crossref
                            source = LightXML.attribute(ref, "dbSource")
                            if source == _name(db_from) &&
                               LightXML.attribute(ref, "dbAccessionId") == id_from
                                key_data = _get_nullable_attribute(ref, "dbResNum")
                            end
                            if source == _name(db_to) &&
                               LightXML.attribute(ref, "dbAccessionId") == id_to
                                value_data = _get_nullable_attribute(ref, "dbResNum")
                            end
                            if !in_chain && source == "PDB" # XML: <crossRefDb dbSource="PDB" ... dbChainId="E"/>
                                in_chain = LightXML.attribute(ref, "dbChainId") == chain
                            end
                        end
                        if !ismissing(key_data) && !ismissing(value_data) && in_chain
                            key = key_data
                            if haskey(mapping, key)
                                @warn string(
                                    "$key is already in the mapping with the value ",
                                    mapping[key],
                                    ". The value is replaced by ",
                                    value_data,
                                )
                            end
                            mapping[key] = value_data
                        end
                    end
                end
            end
        end
    finally
        LightXML.free(xdoc)
    end
    sizehint!(mapping, length(mapping))
end

"""
`parse_file(document::LightXML.XMLDocument, ::Type{SIFTSXML}; chain=All, missings::Bool=true)`

Returns a `Vector{SIFTSResidue}` parsed from a `SIFTSXML` file.
By default, parses all the `chain`s and includes missing residues.
"""
function Utils.parse_file(
    document::LightXML.XMLDocument,
    ::Type{SIFTSXML};
    chain::Union{Type{All},String} = All,
    missings::Bool = true,
)
    vector = SIFTSResidue[]
    for entity in _get_entities(document)
        for segment in _get_segments(entity)
            residues = _get_residues(segment)
            for residue in residues
                missing_residue, sscode, ssname = _get_details(residue)
                if missings || !missing_residue
                    sifts_res = SIFTSResidue(residue, missing_residue, sscode, ssname)
                    if _is_All(chain) ||
                       (!ismissing(sifts_res.PDB) && sifts_res.PDB.chain == chain)
                        push!(vector, sifts_res)
                    end
                end
            end
        end
    end
    vector
end

function Utils.parse_file(fh::Union{IO,AbstractString}, ::Type{SIFTSXML}; kwargs...)
    throw(
        ArgumentError("The SIFTS XML file should have the .xml or the .xml.gz extension."),
    )
end

# Find SIFTSResidue
# -----------------

for F in (:findall, :filter!, :filter)
    @eval begin
        function Base.$(F)(
            f::Function,
            list::AbstractVector{SIFTSResidue},
            db::Type{T},
        ) where {T<:DataBase}
            $(F)(list) do res
                database = get(res, db)
                if !ismissing(database)
                    f(database)
                end
            end
        end
    end
end
