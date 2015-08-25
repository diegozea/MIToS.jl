abstract DataBase{ResNumType}

@auto_hash_equals immutable dbPDBe <: DataBase{Int}
  resnum::Int # Cross referenced residue number
  resname::ASCIIString # Cross referenced residue name
end

@auto_hash_equals immutable dbInterPro <: DataBase{ASCIIString}
  id::ASCIIString
  resnum::ASCIIString # Cross referenced residue number
  resname::ASCIIString # Cross referenced residue name
  evidence::ASCIIString
end

for ref_type in [:dbUniProt, :dbPfam, :dbNCBI]
  @eval begin

    @auto_hash_equals immutable $(ref_type) <: DataBase{Int}
      id::ASCIIString # The cross reference database identifier
      resnum::Int # Cross referenced residue number
      resname::ASCIIString # Cross referenced residue name
      end

  end
end

for ref_type in [:dbPDB, :dbCATH, :dbSCOP]
  @eval begin

    @auto_hash_equals immutable $(ref_type) <: DataBase{ASCIIString}
      id::ASCIIString
      resnum::ASCIIString
      resname::ASCIIString
      chain::ASCIIString
    end

  end
end

name(db::dbPDBe ) = "PDBe"
name(db::dbUniProt) = "UniProt"
name(db::dbPfam) = "Pfam"
name(db::dbNCBI) = "NCBI"
name(db::dbInterPro) = "InterPro"
name(db::dbPDB ) = "PDB"
name(db::dbSCOP) = "SCOP"
name(db::dbCATH) = "CATH"

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

    call(::Type{$(ref_type)}) = $(ref_type)("", "", "", "")

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

    call(::Type{$(ref_type)}) = $(ref_type)("", -9999, "")

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

call(::Type{dbInterPro}) = dbInterPro("", "", "", "")

function call(::Type{dbPDBe}, map::LightXML.XMLElement)
      dbPDBe(
      parse(Int, _get_attribute(map, "dbResNum")),
      _get_attribute(map, "dbResName")
      )
end

call(::Type{dbPDBe}) = dbPDBe(-9999, "")

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

has(res::SIFTSResidue, db::DataBase) = !isnull(getfield(res, symbol(name(db))))

function getdatabase(res::SIFTSResidue, db::DataBase)
  if has(res, db)
    return( get(getfield(res, symbol(name(db)))) )
  end
  nothing
end

function ischain(res::SIFTSResidue, chain::ASCIIString)
  data = getdatabase(res, dbPDB())
  if data !== nothing
    return(data.chain == chain)
  end
  false
end

function has(res::SIFTSResidue, db::DataBase, id::ASCIIString)
 data = getdatabase(res, db)
 data === nothing ? false : return( data.id == id )
end

function has{T}(res::SIFTSResidue, db::DataBase{T}, id::ASCIIString, coord::T)
  data = getdatabase(res, db)
  if data !== nothing && data.id == id
    return(data.resnum == coord)
  end
  false
end

function getcoordinate(res::SIFTSResidue, db::DataBase, id::ASCIIString)
  data = getdatabase(res, db)
  if data !== nothing && data.id == id
    return(data.resnum)
  end
  nothing
end

function getcoordinate(res::SIFTSResidue, db::DataBase, id::ASCIIString, chain::ASCIIString)
  if ischain(res, chain)
    return(getcoordinate(res, db, id))
  end
  nothing
end

# Mapping Functions
# =================

_parse(::Type{Int}, str) = parse(Int, str)
_parse(::Type{ASCIIString}, str) = ascii(str)
@inline _parse(::Type{ASCIIString}, str::ASCIIString) = str

function siftsmapping{F, T}(filename::ASCIIString,
                            db_from::DataBase{F}, id_from::ASCIIString,
                            db_to::DataBase{T}, id_to::ASCIIString; chain::ASCIIString="all", missings::Bool = true)
  mapping = Dict{F, T}()
  for entity in _get_entities(filename)
    segments = _get_segments(entity)
    for segment in segments
      residues = _get_residues(segment)
  		for residue in residues
        in_chain = chain == "all"
        key_data = name(db_from) == "PDBe" ? Nullable(parse(Int, attribute(residue, "dbResNum"))) : Nullable{F}()
        value_data = name(db_to) == "PDBe" ? Nullable(parse(Int, attribute(residue, "dbResNum"))) : Nullable{T}()
	  	  if missings || !_is_missing(residue)
  			  crossref = get_elements_by_tagname(residue, "crossRefDb")
	  			for ref in crossref
            source = attribute(ref, "dbSource")
		  		  if source == name(db_from) && attribute(ref, "dbAccessionId") == id_from
			  		  key_data = Nullable(_parse(F, attribute(ref, "dbResNum")))
				  	end
            if source == name(db_to) && attribute(ref, "dbAccessionId") == id_to
		  			  value_data = Nullable(_parse(T, attribute(ref, "dbResNum")))
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
