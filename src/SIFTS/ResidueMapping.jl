import Base: call

# dbCoordSys Types
# ================

abstract CoordinateSystem

immutable PDBeCoordinates <: CoordinateSystem
  resnum::Int
end

immutable UniProtCoordinates <: CoordinateSystem
  resnum::Int
end

immutable PDBresnumCoordinates <: CoordinateSystem
  resnum::Int
  inscode::Nullable{Char}
end

abstract DataBase{T <: CoordinateSystem}

immutable PDBe <: DataBase{PDBeCoordinates} end

immutable UniProt  <: DataBase{UniProtCoordinates} end
immutable Pfam     <: DataBase{UniProtCoordinates} end
immutable InterPro <: DataBase{UniProtCoordinates} end
immutable NCBI     <: DataBase{UniProtCoordinates} end

immutable PDB  <: DataBase{PDBresnumCoordinates} end
immutable CATH <: DataBase{PDBresnumCoordinates} end
immutable SCOP <: DataBase{PDBresnumCoordinates} end

name(db::UniProt) = "UniProt"
name(db::Pfam) = "Pfam"
name(db::InterPro) = "InterPro"
name(db::NCBI) = "NCBI"
name(db::PDB ) = "PDB "
name(db::CATH) = "CATH"
name(db::SCOP) = "SCOP"

call(::Type{PDBeCoordinates}, str::ASCIIString) = PDBeCoordinates(parse(Int, str))
call(::Type{UniProtCoordinates}, str::ASCIIString) = UniProtCoordinates(parse(Int, str))
function call(::Type{PDBresnumCoordinates}, str::ASCIIString)
  m = match(r"^(\d+)(\D?)$", str)
  if m === nothing
    throw(ErrorException(string(str, " is not a valid PDBresnum.")))
  end
  resnum = parse(Int, m.captures[1])
  ins = m.captures[2]
  PDBresnumCoordinates(resnum, ins != "" ? Nullable{Char}(ins[1]) : Nullable{Char}() )
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

# function siftsmapping{F <: CoordinateSystem, T <: CoordinateSystem}(filename::ASCIIString, chain::ASCIIString, db_from::Type{DataBase{F}}, id_from::ASCIIString, db_to::Type{DataBase{T}}, id_to::ASCIIString, check_observed::Bool = true)
#   if db_from == "PDBe"
#     sifts_PDBe_mapping(filename, chain, db_to, id_to, check_observed)
#   end
#   mapping = Dict{Int, Int}()
# 	segments = _get_segments(_get_entity(_get_entities(filename), chain))
#   for segment in segments
#     residues = _get_residues(segment)
# 		for residue in residues
#       key = -1
#       value = -1
# 		  if !check_observed || !_is_missing(residue)
# 			  crossref = get_elements_by_tagname(residue, "crossRefDb")
# 				for ref in crossref
# 				  if attribute(ref, "dbSource") == db_from && attribute(ref, "dbAccessionId") == id_from
# 					  key = parse(Int, attribute(ref, "dbResNum"))
# 					end
#           if attribute(ref, "dbSource") == db_to && attribute(ref, "dbAccessionId") == id_to
# 					  value = parse(Int, attribute(ref, "dbResNum"))
# 					end
# 				end
#         if key != -1 && value != -1
#           haskey(mapping, key) && warn(string(key, " is already in the mapping with the value ", mapping[key], ". The value is replaced by ", value))
# 					mapping[key] = value
#         end
# 			end
# 		end
#   end
#   sizehint!(mapping, length(mapping))
# end

function siftsmapping(filename::ASCIIString, chain::ASCIIString, db_from::ASCIIString, id_from::ASCIIString, db_to::ASCIIString, id_to::ASCIIString, check_observed::Bool = true)
  if db_from == "PDBe"
    sifts_PDBe_mapping(filename, chain, db_to, id_to, check_observed)
  end
  mapping = Dict{Int, Int}()
	segments = _get_segments(_get_entity(_get_entities(filename), chain))
  for segment in segments
    residues = _get_residues(segment)
		for residue in residues
      key = -1
      value = -1
		  if !check_observed || !_is_missing(residue)
			  crossref = get_elements_by_tagname(residue, "crossRefDb")
				for ref in crossref
				  if attribute(ref, "dbSource") == db_from && attribute(ref, "dbAccessionId") == id_from
					  key = parse(Int, attribute(ref, "dbResNum"))
					end
          if attribute(ref, "dbSource") == db_to && attribute(ref, "dbAccessionId") == id_to
					  value = parse(Int, attribute(ref, "dbResNum"))
					end
				end
        if key != -1 && value != -1
          haskey(mapping, key) && warn(string(key, " is already in the mapping with the value ", mapping[key], ". The value is replaced by ", value))
					mapping[key] = value
        end
			end
		end
  end
  sizehint!(mapping, length(mapping))
end

function siftsPDBemapping(filename::ASCIIString, chain::ASCIIString, db::ASCIIString, id::ASCIIString, check_observed::Bool = true)
  mapping = Dict{Int, Int}()
	segments = _get_segments(_get_entity(_get_entities(filename), chain))
  for segment in segments
    residues = _get_residues(segment)
		for residue in residues
      if attribute(residue, "dbSource") != "PDBe" # XML: <residue dbSource="PDBe" ...
          continue
      end
		  if !check_observed || !_is_missing(residue)
        key = parse(Int, attribute(residue, "dbResNum")) # XML: <residue dbSource="PDBe" dbCoordSys="PDBe" dbResNum="1" ...
			  crossref = get_elements_by_tagname(residue, "crossRefDb")
				for ref in crossref
          if attribute(ref, "dbSource") == db && attribute(ref, "dbAccessionId") == id
            value = attribute(ref, "dbResNum")
            haskey(mapping, key) && warn(string(key, " is already in the mapping with the value ", mapping[key], ". The value is replaced by ", value))
					  mapping[key] = parse(Int, value)
            break
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
