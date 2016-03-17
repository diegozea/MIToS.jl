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

function show(io::IO, res::SIFTSResidue)
  println(io, "SIFTSResidue", res.missing ?  " (missing) " : "")
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
  println(io, "    InterPro: ",  res.InterPro)
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
               missing)
end

call(::Type{SIFTSResidue}, residue::LightXML.XMLElement) =  SIFTSResidue(residue, _is_missing(residue))

# Mapping Functions
# =================

_parse(::Type{Int}, str) = parse(Int, str)
_parse(::Type{ASCIIString}, str) = ascii(str)
@inline _parse(::Type{ASCIIString}, str::ASCIIString) = str

"""Parses a SIFTS XML file and returns a `Dict` between residue numbers of two `DataBase`s  with the given identifiers.
A `chain` could be specified ("all" by default).
If `missings` is `true` (default) all the residues are used, even if they haven’t coordinates in the PDB file."""
function siftsmapping{F, T}(filename::AbstractString,
                            db_from::Type{F}, id_from::ASCIIString,
                            db_to::Type{T}, id_to::ASCIIString; chain::ASCIIString="all", missings::Bool = true)
  mapping = Dict{_number_type(F), _number_type(T)}()
  xdoc = parse_file(filename)
  try
    for entity in _get_entities(xdoc)
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
  finally
    free(xdoc)
  end
  sizehint!(mapping, length(mapping))
end

function parse(document::LightXML.XMLDocument, ::Type{SIFTSXML}; chain::ASCIIString="all", missings::Bool = true)
  vector = SIFTSResidue[]
  for entity in _get_entities(document)
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

# Find SIFTSResidue
# -----------------

"Returns `true` if the tests are successfully passed for that `DataBase` sub-type on that `SIFTSResidue`."
function isobject{T <: DataBase}(res::SIFTSResidue, ::Type{T}, tests::AbstractTest...)
  dbfield = getfield(res, symbol(name(T)))
  if !isnull(dbfield)
    return(isobject(get(dbfield), tests...))
  end
  false
end

function findobjects{T <: DataBase}(vector::AbstractVector{SIFTSResidue}, ::Type{T}, tests::AbstractTest...)
  dbname = symbol(name(T))
  len = length(vector)
  indexes = Array(Int, len)
  n = 0
  for i in 1:len
    dbfield = getfield(vector[i], dbname)
    if !isnull(dbfield) && isobject(get(dbfield), tests...)
      n += 1
      indexes[n] = i
    end
  end
  resize!(indexes, n)
end

function collectobjects{T <: DataBase}(vector::AbstractVector{SIFTSResidue}, ::Type{T}, tests::AbstractTest...)
  dbname = symbol(name(T))
  len = length(vector)
  elements = Array(SIFTSResidue, len)
  n = 0
  for i in 1:len
    dbfield = getfield(vector[i], dbname)
    if !isnull(dbfield) && isobject(get(dbfield), tests...)
      n += 1
      elements[n] = vector[i]
    end
  end
  resize!(elements, n)
end

"""
`capture(res::SIFTSResidue, db_capture, field, db_test, tests...)`

Takes a `SIFTSResidue`, a `db...` type and a `Symbol` with the name of the `field` to capture from that database.
Returns a `Nullable` with the field content if the `tests` are passed over a determined database (`db_test`).
The function Returns `nothing` if the DataBase to test or capture is null in the `SIFTSResidue`.
"""
function capture{C <: DataBase, T <: DataBase}(res::SIFTSResidue, db_capture::Type{C}, field::Symbol, db_test::Type{T}, tests::AbstractTest...)
  dbfield_capture = getfield(res, symbol(name(C)))
  dbfield_test = getfield(res, symbol(name(T)))
  if !isnull(dbfield_capture) && !isnull(dbfield_test)
    captured = getfield(get(dbfield_capture), field)
    for test in tests
      if !Utils._test(get(dbfield_test), test)
        return(Nullable{typeof(captured)}())
      end
    end
    return(Nullable(captured))
  end
  nothing
end

"Returns the **type** of the `field` in the first object of the collection"
function guess_type{T <: DataBase}(collection::AbstractVector{SIFTSResidue}, ::Type{T}, field::Symbol)
  for res in collection
    dbfield = getfield(res, symbol(name(T)))
    if !isnull(dbfield)
      data = get(dbfield)
      if field in fieldnames(data)
        return(fieldtype(typeof(data), field))
      end
    end
  end
  Any
end

"""
`collectcaptures(vector::AbstractVector{SIFTSResidue}, db_capture, field, db_test, tests)`

Calls the `capture` function and returns a vector of `Nullable`s with the captured `field`s.
The element is null if any test fails or the object hasn't the `field`.
"""
function collectcaptures{C <: DataBase, T <: DataBase}(vector::AbstractVector{SIFTSResidue}, db_capture::Type{C}, field::Symbol, db_test::Type{T}, tests::AbstractTest...)
  Nullable{guess_type(vector, C, field)}[ capture(object, C, field, T, tests...) for object in vector ]
end
