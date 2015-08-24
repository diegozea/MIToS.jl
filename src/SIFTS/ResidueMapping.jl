import Base: call, string, print, write, show, convert, isnull

# dbCoordSys Types
# ================

abstract CoordinateSystem

@auto_hash_equals immutable PDBeCoordinate <: CoordinateSystem
  resnum::Int
end

@auto_hash_equals immutable UniProtCoordinate <: CoordinateSystem
  resnum::Int
end

@auto_hash_equals immutable PDBresnumCoordinate <: CoordinateSystem
  resnum::Int
  inscode::ASCIIString
end

call(::Type{PDBresnumCoordinate}, resnum::Int) = PDBresnumCoordinate(resnum, "")
call(::Type{PDBresnumCoordinate}, resnum::Int, ins::Char) = PDBresnumCoordinate(resnum, string(ins))

convert(::Type{PDBeCoordinate}, str::ASCIIString) = PDBeCoordinate(parse(Int, str))
convert(::Type{UniProtCoordinate}, str::ASCIIString) = UniProtCoordinate(parse(Int, str))
function convert(::Type{PDBresnumCoordinate}, str::ASCIIString)
  m = match(r"^(-?\d+)(\D?)$", str)
  if m === nothing
    throw(ErrorException(string(str, " is not a valid PDBresnum.")))
  end
  PDBresnumCoordinate(parse(Int, m.captures[1]), m.captures[2])
end

call(::Type{PDBeCoordinate}, str::ASCIIString) = PDBeCoordinate(parse(Int, str))
call(::Type{UniProtCoordinate}, str::ASCIIString) = UniProtCoordinate(parse(Int, str))
call(::Type{PDBresnumCoordinate}, str::ASCIIString) = convert(PDBresnumCoordinate, str)

for (typ, sys) in [(:PDBeCoordinate, "PDBe"), (:UniProtCoordinate, "UniProt"), (:PDBresnumCoordinate, "PDBresnum")]
  @eval begin
    function call(::Type{$(typ)}, elem::LightXML.XMLElement)
      system = attribute(elem, "dbCoordSys")
      if system == $sys
        return(convert($(typ), attribute(elem, "dbResNum")))
      else
        throw(ErrorException(string($sys, " is the expected co-ordinate system, but the system is ", system)))
      end
    end
  end
end

string(coord::Union{PDBeCoordinate, UniProtCoordinate}) = string(coord.resnum)
string(coord::PDBresnumCoordinate) = string(coord.resnum, coord.inscode)

convert(::Type{Int}, coord::CoordinateSystem) = coord.resnum
isnull(coord::CoordinateSystem) = Int(coord) == -9999

show{T <: CoordinateSystem}(io::IO, coord::T) = print(io, string(T, "(\"", coord, "\")"))
print(io::IO, coord::CoordinateSystem) = print(io, string(coord))
write(io::IO, coord::CoordinateSystem) = write(io, string(coord))

abstract DataBase{T <: CoordinateSystem}

@auto_hash_equals immutable RefPDBe <: DataBase{PDBeCoordinate}
  dbResNum::PDBeCoordinate # Cross referenced residue number
  dbResName::ASCIIString # Cross referenced residue name
end

# @auto_hash_equals immutable RefInterPro <: DataBase{UniProtCoordinate}
#   dbAccessionId::ASCIIString
#   dbResNum::UniProtCoordinate # Cross referenced residue number
#   dbResName::ASCIIString # Cross referenced residue name
#   dbEvidence::ASCIIString
# end

for ref_type in [:RefUniProt, :RefPfam, :RefNCBI]
  @eval begin

    @auto_hash_equals immutable $(ref_type) <: DataBase{UniProtCoordinate}
      dbAccessionId::ASCIIString # The cross reference database identifier
      dbResNum::UniProtCoordinate # Cross referenced residue number
      dbResName::ASCIIString # Cross referenced residue name
      end

  end
end

for ref_type in [:RefPDB, :RefCATH, :RefSCOP]
  @eval begin

    @auto_hash_equals immutable $(ref_type) <: DataBase{PDBresnumCoordinate}
      dbAccessionId::ASCIIString
      dbResNum::PDBresnumCoordinate
      dbResName::ASCIIString
      dbChainId::ASCIIString
    end

  end
end

const dbPDBe = RefPDBe("-9999", "")

const dbUniProt = RefUniProt("", "-9999", "")
const dbPfam = RefPfam("", "-9999", "")
const dbNCBI = RefNCBI("", "-9999", "")
#const dbInterPro = RefInterPro("", "-9999", "", "")

const dbPDB = RefPDB("", "-9999", "", "")
const dbSCOP = RefSCOP("", "-9999", "", "")
const dbCATH = RefCATH("", "-9999", "", "")

name(db::RefPDBe ) = "PDBe"

name(db::RefUniProt) = "UniProt"
name(db::RefPfam) = "Pfam"
name(db::RefNCBI) = "NCBI"
#name(db::RefInterPro) = "InterPro"

name(db::RefPDB ) = "PDB"
name(db::RefSCOP) = "SCOP"
name(db::RefCATH) = "CATH"

"""Returns "" if the attributte is missing"""
function _get_attribute(elem::LightXML.XMLElement, attr::ASCIIString)
  text = attribute(elem, attr)
  if text === nothing
    return("")
  else
    return(text)
  end
end

for ref_type in [:RefPDB, :RefCATH, :RefSCOP]
  @eval begin
    function call(::Type{$(ref_type)}, map::LightXML.XMLElement)
      $(ref_type)(
      _get_attribute(map, "dbAccessionId"),
      PDBresnumCoordinate(map),
      _get_attribute(map, "dbResName"),
      _get_attribute(map, "dbChainId")
      )
    end
  end
end

for ref_type in [:RefUniProt, :RefPfam, :RefNCBI]
  @eval begin
    function call(::Type{$(ref_type)}, map::LightXML.XMLElement)
      $(ref_type)(
      _get_attribute(map, "dbAccessionId"),
      UniProtCoordinate(map),
      _get_attribute(map, "dbResName")
      )
    end
  end
end

# function call(::Type{RefInterPro}, map::LightXML.XMLElement)
#   if _get_attribute(map, "dbResNum") == "None"
#     warn("There is not dbResNum for InterPro, using -9999")
#     resnum = UniProtCoordinate("-9999")
#   else
#     resnum = UniProtCoordinate(map)
#   end
#   RefInterPro(
#     _get_attribute(map, "dbAccessionId"),
#     resnum,
#     _get_attribute(map, "dbResName"),
#     _get_attribute(map, "dbEvidence")
#     )
# end

function call(::Type{RefPDBe}, map::LightXML.XMLElement)
      RefPDBe(
      PDBeCoordinate(map),
      _get_attribute(map, "dbResName")
      )
end

@auto_hash_equals immutable SIFTSResidue
  PDBe::RefPDBe
  UniProt::Nullable{RefUniProt}
  Pfam::Nullable{RefPfam}
  NCBI::Nullable{RefNCBI}
#  InterPro::Array{RefInterPro,1}
  PDB::Nullable{RefPDB}
  SCOP::Nullable{RefSCOP}
  CATH::Nullable{RefCATH}
  missing::Bool # XML: <residueDetail dbSource="PDBe" property="Annotation" ...
end

function call(::Type{SIFTSResidue}, residue::LightXML.XMLElement, missing::Bool)
  PDBe = RefPDBe(residue)
  UniProt = Nullable{RefUniProt}()
  Pfam = Nullable{RefPfam}()
  NCBI = Nullable{RefNCBI}()
#  InterPro = RefInterPro[]
  PDB = Nullable{RefPDB}()
  SCOP = Nullable{RefSCOP}()
  CATH = Nullable{RefCATH}()
  for crossref in get_elements_by_tagname(residue, "crossRefDb")
    db = attribute(crossref, "dbSource")
    if db == "UniProt"
      UniProt = Nullable(RefUniProt(crossref))
    elseif db == "Pfam"
      Pfam = Nullable(RefPfam(crossref))
    elseif db == "NCBI"
      NCBI = Nullable(RefNCBI(crossref))
#    elseif db == "InterPro"
#      push!(InterPro, RefInterPro(crossref))
    elseif db == "PDB"
      PDB = Nullable(RefPDB(crossref))
    elseif db == "SCOP"
      SCOP = Nullable(RefSCOP(crossref))
    elseif db == "CATH"
      CATH = Nullable(RefCATH(crossref))
    else
      # InterPro isn't supported rigth now, ResNum can be "None" instead of Int
      db != "InterPro" && warn(string(db, " is not in the MIToS' DataBases."))
    end
  end
  SIFTSResidue(PDBe,
               UniProt,
               Pfam,
               NCBI,
#               InterPro,
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
  data = getdatabase(res, dbPDB)
  if data !== nothing
    return(data.dbChainId == chain)
  end
  false
end

function has(res::SIFTSResidue, db::DataBase, id::ASCIIString)
 data = getdatabase(res, db)
 data === nothing ? false : return( data.dbAccessionId == id )
end

function has{T}(res::SIFTSResidue, db::DataBase{T}, id::ASCIIString, coord::T)
  data = getdatabase(res, db)
  if data !== nothing && data.dbAccessionId == id
    return(data.dbResNum == coord)
  end
  false
end

function getcoordinate{T}(res::SIFTSResidue, db::DataBase{T}, id::ASCIIString)
  data = getdatabase(res, db)
  if data !== nothing && data.dbAccessionId == id
    return(data.dbResNum)
  end
  nothing
end

function getcoordinate{T}(res::SIFTSResidue, db::DataBase{T}, id::ASCIIString, chain::ASCIIString)
  if ischain(res, chain)
    return(getcoordinate(res, db, id))
  end
  nothing
end

# Mapping Functions
# =================

function siftsmapping{F <: CoordinateSystem, T <: CoordinateSystem}(filename::ASCIIString,
                                                              db_from::DataBase{F}, id_from::ASCIIString,
                                                              db_to::DataBase{T}, id_to::ASCIIString; chain::ASCIIString="all", missings::Bool = true)
  mapping = Dict{F, T}()
  for entity in _get_entities(filename)
    segments = _get_segments(entity)
    for segment in segments
      residues = _get_residues(segment)
  		for residue in residues
        in_chain = chain == "all"
        key_data = name(db_from) == "PDBe" ? Nullable(PDBeCoordinate(residue)) : Nullable{F}()
        value_data = name(db_to) == "PDBe" ? Nullable(PDBeCoordinate(residue)) : Nullable{F}()
	  	  if missings || !_is_missing(residue)
  			  crossref = get_elements_by_tagname(residue, "crossRefDb")
	  			for ref in crossref
            source = attribute(ref, "dbSource")
		  		  if source == name(db_from) && attribute(ref, "dbAccessionId") == id_from
			  		  key_data = Nullable(F(ref))
				  	end
            if source == name(db_to) && attribute(ref, "dbAccessionId") == id_to
		  			  value_data = Nullable(T(ref))
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
          if chain == "all" || (!isnull(sifts_res.PDB) && get(sifts_res.PDB).dbChainId == chain)
            push!(vector, sifts_res)
          end
			  end
		  end
    end
  end
  vector
end
