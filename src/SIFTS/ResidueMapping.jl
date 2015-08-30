abstract DataBase

@auto_hash_equals immutable dbPDBe <: DataBase
  number::Int # Cross referenced residue number
  name::ASCIIString # Cross referenced residue name
end

@inline _number_type(::Type{dbPDBe}) = Int

@auto_hash_equals immutable dbInterPro <: DataBase
  id::ASCIIString
  number::ASCIIString # Cross referenced residue number
  name::ASCIIString # Cross referenced residue name
  evidence::ASCIIString
end

@inline _number_type(::Type{dbInterPro}) = ASCIIString

for ref_type in [:dbUniProt, :dbPfam, :dbNCBI]
  @eval begin

    @auto_hash_equals immutable $(ref_type) <: DataBase
      id::ASCIIString # The cross reference database identifier
      number::Int # Cross referenced residue number
      name::ASCIIString # Cross referenced residue name
      end

    @inline _number_type(::Type{$(ref_type)}) = Int

  end
end

for ref_type in [:dbPDB, :dbCATH, :dbSCOP]
  @eval begin

    @auto_hash_equals immutable $(ref_type) <: DataBase
      id::ASCIIString
      number::ASCIIString
      name::ASCIIString
      chain::ASCIIString
    end

    @inline _number_type(::Type{$(ref_type)}) = ASCIIString

  end
end

@inline name(::Type{dbPDBe})    = "PDBe"
@inline name(::Type{dbUniProt}) = "UniProt"
@inline name(::Type{dbPfam})    = "Pfam"
@inline name(::Type{dbNCBI})    = "NCBI"
@inline name(::Type{dbInterPro})= "InterPro"
@inline name(::Type{dbPDB})     = "PDB"
@inline name(::Type{dbSCOP})    = "SCOP"
@inline name(::Type{dbCATH})    = "CATH"

"""Returns "" if the attributte is missing"""
function _get_attribute(elem::LightXML.XMLElement, attr::ASCIIString)
  text = attribute(elem, attr)
  if text === nothing
    return("")
  else
    return(text)
  end
end

for ref_type in [:dbPDB, :dbCATH, :dbSCOP]
  @eval begin

    function call(::Type{$(ref_type)}, map::LightXML.XMLElement)
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

    function call(::Type{$(ref_type)}, map::LightXML.XMLElement)
      $(ref_type)(
      _get_attribute(map, "dbAccessionId"),
      parse(Int, _get_attribute(map, "dbResNum")),
      _get_attribute(map, "dbResName")
      )
    end

  end
end

function call(::Type{dbInterPro}, map::LightXML.XMLElement)
  dbInterPro(
    _get_attribute(map, "dbAccessionId"),
    _get_attribute(map, "dbResNum"),
    _get_attribute(map, "dbResName"),
    _get_attribute(map, "dbEvidence")
    )
end

function call(::Type{dbPDBe}, map::LightXML.XMLElement)
      dbPDBe(
      parse(Int, _get_attribute(map, "dbResNum")),
      _get_attribute(map, "dbResName")
      )
end

@auto_hash_equals immutable SIFTSResidue
  PDBe::dbPDBe
  UniProt::Nullable{dbUniProt}
  Pfam::Nullable{dbPfam}
  NCBI::Nullable{dbNCBI}
  InterPro::Array{dbInterPro, 1}
  PDB::Nullable{dbPDB}
  SCOP::Nullable{dbSCOP}
  CATH::Nullable{dbCATH}
  missing::Bool # XML: <residueDetail dbSource="PDBe" property="Annotation" ...
end

function call(::Type{SIFTSResidue}, residue::LightXML.XMLElement, missing::Bool)
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
      db != "InterPro" && warn(string(db, " is not in the MIToS' DataBases."))
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
               missing)
end

call(::Type{SIFTSResidue}, residue::LightXML.XMLElement) =  SIFTSResidue(residue, _is_missing(residue))

# Asking to SIFTSResidue
# ======================

has{T <: DataBase}(res::SIFTSResidue, ::Type{T}) = !isnull(getfield(res, symbol(name(T))))

function getdatabase{T <: DataBase}(res::SIFTSResidue, ::Type{T})
  if has(res, T)
    return( get(getfield(res, symbol(name(T)))) )
  end
  nothing
end

function ischain(res::SIFTSResidue, chain::ASCIIString)
  data = getdatabase(res, dbPDB)
  if data !== nothing
    return(data.chain == chain)
  end
  false
end

function has{T <: DataBase}(res::SIFTSResidue, ::Type{T}, id::ASCIIString)
 data = getdatabase(res, T)
 data === nothing ? false : return( data.id == id )
end

function has{T <: DataBase}(res::SIFTSResidue, ::Type{T}, id::ASCIIString, coord)
  data = getdatabase(res, T)
  if data !== nothing && data.id == id
    return(data.number == coord)
  end
  false
end

function getcoordinate{T <: DataBase}(res::SIFTSResidue, ::Type{T}, id::ASCIIString)
  data = getdatabase(res, T)
  if data !== nothing && data.id == id
    return(data.number)
  end
  nothing
end

function getcoordinate{T <: DataBase}(res::SIFTSResidue, ::Type{T}, id::ASCIIString, chain::ASCIIString)
  if ischain(res, chain)
    return(getcoordinate(res, T, id))
  end
  nothing
end

# Mapping Functions
# =================

_parse(::Type{Int}, str) = parse(Int, str)
_parse(::Type{ASCIIString}, str) = ascii(str)
@inline _parse(::Type{ASCIIString}, str::ASCIIString) = str

function siftsmapping{F, T}(filename::ASCIIString,
                            db_from::Type{F}, id_from::ASCIIString,
                            db_to::Type{T}, id_to::ASCIIString; chain::ASCIIString="all", missings::Bool = true)
  mapping = Dict{_number_type(F), _number_type(T)}()
  for entity in _get_entities(filename)
    segments = _get_segments(entity)
    for segment in segments
      residues = _get_residues(segment)
  		for residue in residues
        in_chain = chain == "all"
        key_data = name(db_from) == "PDBe" ? Nullable(parse(Int, attribute(residue, "dbResNum"))) : Nullable{_number_type(F)}()
        value_data = name(db_to) == "PDBe" ? Nullable(parse(Int, attribute(residue, "dbResNum"))) : Nullable{_number_type(T)}()
	  	  if missings || !_is_missing(residue)
  			  crossref = get_elements_by_tagname(residue, "crossRefDb")
	  			for ref in crossref
            source = attribute(ref, "dbSource")
		  		  if source == name(db_from) && attribute(ref, "dbAccessionId") == id_from
			  		  key_data = Nullable(_parse(_number_type(F), attribute(ref, "dbResNum")))
				  	end
            if source == name(db_to) && attribute(ref, "dbAccessionId") == id_to
		  			  value_data = Nullable(_parse(_number_type(T), attribute(ref, "dbResNum")))
			  		end
            if !in_chain && source == "PDB" # XML: <crossRefDb dbSource="PDB" ... dbChainId="E"/>
              in_chain = attribute(ref, "dbChainId") == chain
            end
		  		end
          if !isnull(key_data) && !isnull(value_data) && in_chain
            key = get(key_data)
            haskey(mapping, key) && warn(string(key, " is already in the mapping with the value ", mapping[key], ". The value is replaced by ", get(value_data)))
			  		mapping[key] = get(value_data)
          end
			  end
		  end
    end
  end
  sizehint!(mapping, length(mapping))
end

function siftsresidues(filename::ASCIIString; chain::ASCIIString="all", missings::Bool = true)
  vector = SIFTSResidue[]
  for entity in _get_entities(filename)
    for segment in _get_segments(entity)
      residues = _get_residues(segment)
		  for residue in residues
        missing = _is_missing(residue)
		    if missings || !missing
          sifts_res = SIFTSResidue(residue, missing)
          if chain == "all" || (!isnull(sifts_res.PDB) && get(sifts_res.PDB).chain == chain)
            push!(vector, sifts_res)
          end
			  end
		  end
    end
  end
  vector
end
