import Base: call, string, print, write, show, convert, hash, ==

# dbCoordSys Types
# ================

abstract CoordinateSystem

immutable PDBeCoordinate <: CoordinateSystem
  resnum::Int
end

immutable UniProtCoordinate <: CoordinateSystem
  resnum::Int
end

immutable PDBresnumCoordinate <: CoordinateSystem
  resnum::Int
  inscode::Nullable{Char}
end

string(coord::Union{PDBeCoordinate, UniProtCoordinate}) = string(coord.resnum)
function string(coord::PDBresnumCoordinate)
  if isnull(coord.inscode)
    return(string(coord.resnum))
  else
    return(string(coord.resnum, get(coord.inscode)))
  end
end

convert(::Type{Int}, coord::CoordinateSystem) = coord.resnum

show{T <: CoordinateSystem}(io::IO, coord::T) = print(io, string(T, "(\"", coord, "\")"))
print(io::IO, coord::CoordinateSystem) = print(io, string(coord))
write(io::IO, coord::CoordinateSystem) = write(io, string(coord))

abstract DataBase{T <: CoordinateSystem}

immutable dbPDBe <: DataBase{PDBeCoordinate} end

immutable dbUniProt  <: DataBase{UniProtCoordinate} end
immutable dbPfam     <: DataBase{UniProtCoordinate} end
immutable dbInterPro <: DataBase{UniProtCoordinate} end
immutable dbNCBI     <: DataBase{UniProtCoordinate} end

immutable dbPDB  <: DataBase{PDBresnumCoordinate} end
immutable dbCATH <: DataBase{PDBresnumCoordinate} end
immutable dbSCOP <: DataBase{PDBresnumCoordinate} end

name(db::dbPDBe ) = "PDBe"

name(db::dbUniProt) = "UniProt"
name(db::dbPfam) = "Pfam"
name(db::dbInterPro) = "InterPro"
name(db::dbNCBI) = "NCBI"

name(db::dbPDB ) = "PDB"
name(db::dbCATH) = "CATH"
name(db::dbSCOP) = "SCOP"

call(::Type{PDBeCoordinate}, str::ASCIIString) = PDBeCoordinate(parse(Int, str))
call(::Type{UniProtCoordinate}, str::ASCIIString) = UniProtCoordinate(parse(Int, str))
call(::Type{PDBresnumCoordinate}, resnum::Int) = PDBresnumCoordinate(resnum, Nullable{Char}())
call(::Type{PDBresnumCoordinate}, resnum::Int, ins::AbstractString) = PDBresnumCoordinate(resnum, ins != "" ? Nullable{Char}(ins[1]) : Nullable{Char}())
function call(::Type{PDBresnumCoordinate}, str::ASCIIString)
  m = match(r"^(-?\d+)(\D?)$", str)
  if m === nothing
    throw(ErrorException(string(str, " is not a valid PDBresnum.")))
  end
  PDBresnumCoordinate(parse(Int, m.captures[1]), m.captures[2])
end

function ==(x::PDBresnumCoordinate, y::PDBresnumCoordinate)
  if x.resnum == y.resnum
    x_ins_flag = isnull(x.inscode)
    y_ins_flag = isnull(y.inscode)
    if x_ins_flag && y_ins_flag
      return(true)
    elseif !x_ins_flag && !y_ins_flag
      return(get(x.inscode) == get(x.inscode))
    end
  end
  false
end

hash(x::CoordinateSystem, h) = hash(string(x), h)

for (typ, sys) in [(:PDBeCoordinate, "PDBe"), (:UniProtCoordinate, "UniProt"), (:PDBresnumCoordinate, "PDBresnum")]
  @eval begin
    function call(::Type{$(typ)}, elem::LightXML.XMLElement)
      system = attribute(elem, "dbCoordSys")
      if system == $sys
        return($(typ)(attribute(elem, "dbResNum")))
      else
        throw(ErrorException(string($sys, " is the expected co-ordinate system, but the system is ", system)))
      end
    end
  end
end

# Mapping Types
# =============

immutable crossRefDb
  dbSource::ASCIIString # The database being cross referenced
  dbVersion::ASCIIString # The cross referenced database version
  dbCoordSys::ASCIIString # The cross referenced database co-ordinate system
  dbAccessionId::ASCIIString # The cross reference database identifier
  dbResNum::ASCIIString # Cross referenced residue number
  dbResName::ASCIIString # Cross referenced residue name
  dbChainId # Cross referenced chain id
end

"""Returns "" if the attributte is missing"""
function _get_attribute(elem::LightXML.XMLElement, attr::ASCIIString)
  text = attribute(elem, attr)
  if text === nothing
    return("")
  else
    return(text)
  end
end

function call(::Type{crossRefDb}, map::LightXML.XMLElement)
  crossRefDb(
    _get_attribute(map, "dbSource"),
    _get_attribute(map, "dbVersion"),
    _get_attribute(map, "dbCoordSys"),
    _get_attribute(map, "dbAccessionId"),
    _get_attribute(map, "dbResNum"),
    _get_attribute(map, "dbResName"),
    _get_attribute(map, "dbChainId")
    )
end

immutable SIFTSResidue
  db::ASCIIString # dbSource
  number::ASCIIString # dbResNum: The residue number
  name::ASCIIString # dbResName: The residue name
  mapping::Array{crossRefDb,1} # XML: <crossRefDb ...
  missing::Bool # XML: <residueDetail dbSource="PDBe" property="Annotation" ...
end

function call(::Type{SIFTSResidue}, residue::LightXML.XMLElement)
  SIFTSResidue(
    _get_attribute(residue, "dbSource"),
    _get_attribute(residue, "dbResNum"),
    _get_attribute(residue, "dbResName"),
    crossRefDb[ crossRefDb(map) for map in get_elements_by_tagname(residue, "crossRefDb")],
    _is_missing(residue)
  )
end

function call(::Type{SIFTSResidue}, residue::LightXML.XMLElement, missing::Bool)
  SIFTSResidue(
    _get_attribute(residue, "dbSource"),
    _get_attribute(residue, "dbResNum"),
    _get_attribute(residue, "dbResName"),
    crossRefDb[ crossRefDb(map) for map in get_elements_by_tagname(residue, "crossRefDb")],
    missing
  )
end

# Mapping Functions
# =================

function siftsmapping{F <: CoordinateSystem, T <: CoordinateSystem}(filename::ASCIIString, chain::ASCIIString,
                                                              db_from::DataBase{F}, id_from::ASCIIString,
                                                              db_to::DataBase{T}, id_to::ASCIIString, check_observed::Bool = true)
  mapping = Dict{F, T}()
	segments = _get_segments(_get_entity(_get_entities(filename), chain))
  for segment in segments
    residues = _get_residues(segment)
		for residue in residues
      key_data = Nullable{F}
      value_data = Nullable{T}
		  if !check_observed || !_is_missing(residue)
			  crossref = get_elements_by_tagname(residue, "crossRefDb")
				for ref in crossref
				  if attribute(ref, "dbSource") == name(db_from) && attribute(ref, "dbAccessionId") == id_from
					  key_data = Nullable(F(ref))
					end
          if attribute(ref, "dbSource") == name(db_to) && attribute(ref, "dbAccessionId") == id_to
					  value_data = Nullable(T(ref))
					end
				end
        if !isnull(key_data) && !isnull(value_data)
          key = get(key_data)
          haskey(mapping, key) && warn(string(key, " is already in the mapping with the value ", mapping[key], ". The value is replaced by ", get(value_data)))
					mapping[key] = get(value_data)
        end
			end
		end
  end
  sizehint!(mapping, length(mapping))
end

immutable dbPDBe <: DataBase{PDBeCoordinate} end

function siftsmapping{T <: CoordinateSystem}(filename::ASCIIString, chain::ASCIIString, db_from::dbPDBe,
                                             db_to::DataBase{T}, id_to::ASCIIString, check_observed::Bool = true)
  mapping = Dict{PDBeCoordinate, T}()
	segments = _get_segments(_get_entity(_get_entities(filename), chain))
  for segment in segments
    residues = _get_residues(segment)
		for residue in residues
      if attribute(residue, "dbSource") != "PDBe" # XML: <residue dbSource="PDBe" ...
          continue
      end
      key = PDBeCoordinate(residue)
		  if !check_observed || !_is_missing(residue)
			  crossref = get_elements_by_tagname(residue, "crossRefDb")
				for ref in crossref
				  if attribute(ref, "dbSource") == name(db_to) && attribute(ref, "dbAccessionId") == id_to
            value = T(ref)
            haskey(mapping, key) && warn(string(key, " is already in the mapping with the value ", mapping[key], ". The value is replaced by ", value))
					  mapping[key] = value
					end
				end
			end
		end
  end
  sizehint!(mapping, length(mapping))
end

function siftsmapping{F <: CoordinateSystem}(filename::ASCIIString, chain::ASCIIString, db_from::DataBase{F}, id_from::ASCIIString,
                                             db_to::dbPDBe, check_observed::Bool = true)
  mapping = Dict{F, PDBeCoordinate}()
	segments = _get_segments(_get_entity(_get_entities(filename), chain))
  for segment in segments
    residues = _get_residues(segment)
		for residue in residues
      if attribute(residue, "dbSource") != "PDBe" # XML: <residue dbSource="PDBe" ...
          continue
      end
      value = PDBeCoordinate(residue)
		  if !check_observed || !_is_missing(residue)
			  crossref = get_elements_by_tagname(residue, "crossRefDb")
				for ref in crossref
				  if attribute(ref, "dbSource") == name(db_from) && attribute(ref, "dbAccessionId") == id_from
            key = F(ref, "dbResNum")
            haskey(mapping, key) && warn(string(key, " is already in the mapping with the value ", mapping[key], ". The value is replaced by ", value))
					  mapping[key] = value
					end
				end
			end
		end
  end
  sizehint!(mapping, length(mapping))
end

function siftsresidues(filename::ASCIIString, chain::ASCIIString, check_observed::Bool = true)
  vector = SIFTSResidue[]
	segments = _get_segments(_get_entity(_get_entities(filename), chain))
  for segment in segments
    residues = _get_residues(segment)
		for residue in residues
      missing = _is_missing(residue)
		  if !check_observed || !missing
        push!(vector, SIFTSResidue(residue, missing))
			end
		end
  end
  vector
end
