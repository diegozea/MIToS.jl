import Base: call

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
  number::ASCIIString # dbResNum: The residue number
  name::ASCIIString # dbResName: The residue name
  mapping::Array{crossRefDb,1} # XML: <crossRefDb ...
  missing::Bool # XML: <residueDetail dbSource="PDBe" property="Annotation" ...
end

function call(::Type{SIFTSResidue}, residue::LightXML.XMLElement)
  SIFTSResidue(
    _get_attribute(residue, "dbResNum"),
    _get_attribute(residue, "dbResName"),
    crossRefDb[ crossRefDb(map) for map in get_elements_by_tagname(residue, "crossRefDb")],
    _is_missing(residue)
  )
end

function call(::Type{SIFTSResidue}, residue::LightXML.XMLElement, missing::Bool)
  SIFTSResidue(
    _get_attribute(residue, "dbResNum"),
    _get_attribute(residue, "dbResName"),
    crossRefDb[ crossRefDb(map) for map in get_elements_by_tagname(residue, "crossRefDb")],
    missing
  )
end

# Mapping Functions
# =================

function siftsmapping(filename::ASCIIString, chain::ASCIIString, dbid::ASCIIString; db::ASCIIString = "Pfam", check_observed::Bool = true)
  mapping = Dict{Int,Int}()
	segments = _get_segments(_get_entity(_get_entities(filename), chain))
  for segment in segments
    residues = _get_residues(segment)
		for residue in residues
		  if !check_observed || !_is_missing(residue)
        key = parse(Int, attribute(residue, "dbResNum")) # XML: <residue dbSource="PDBe" ...
			  crossref = get_elements_by_tagname(residue, "crossRefDb")
				for ref in crossref
				  if attribute(ref, "dbSource") == db && attribute(ref, "dbAccessionId") == dbid
					  mapping[key] = parse(Int, attribute(ref, "dbResNum"))
            break
					end
				end
			end
		end
  end
  sizehint!(mapping, length(mapping))
end

function siftsresidues(filename::ASCIIString, chain::ASCIIString; check_observed::Bool = true)
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
