# Download SIFTS
# ==============
"""
Download the gzipped SIFTS xml  for the `pdbcode`.
The extension of the downloaded file is `.xml.gz` by default. The `filename` can be changed, but the `.gz` at the end is mandatory.
"""
function downloadsifts(pdbcode::ASCIIString; filename::ASCIIString="$(lowercase(pdbcode)).xml.gz")
  if ismatch(r"^\w{4}$"i, pdbcode)
    download(string("ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/split_xml/", lowercase(pdbcode[2:3]), "/", lowercase(pdbcode), ".xml.gz"), filename)
  else
    throw(string(pdbcode, " is not a correct PDB"))
  end
end

# Internal Parser Functions
# =========================

"""
Gets the entities of a SIFTS XML. Each entity is a PDB chain. For example:
```
<entry dbSource="PDBe" ...
  ...
  <entity type="protein" entityId="A">
    ...
  </entity>
  <entity type="protein" entityId="B">
    ...
  </entity>
<\entry>
```
"""
function _get_entities(filename)
  sifts = parse_file(filename)
	siftsroot = root(sifts)
	get_elements_by_tagname(siftsroot, "entity")
end

# """
# Returns the protein entity with PDB `chain` as `entityId` from the array of entities.
# WARNING: Sometimes there are more chains than entities!
# """
# function _get_entity(entities, chain::ASCIIString)
#   for entity in entities
#     # XML: <entity type="protein" entityId="A">
#     if attribute(entity, "entityId") == chain && attribute(entity, "type") == "protein"
#       return(entity)
#     end
#   end
# end

"""
Gets an array of the segments, the continuous region of an entity.
Chimeras and expression tags generates more than one segment for example.
"""
_get_segments(entity) = get_elements_by_tagname(entity, "segment")

"""
Returns the difference between start and end, for example `6` for:
```
<segment segId="2vqc_A_1_7" start="1" end="7">
```
"""
_length_segment(segment) = parse(Int, attribute(segment, "end")) - parse(Int, attribute(segment, "start"))

"""
Returns an `Int` with the index of the largest segment.
For example returns 2 for 2VQC:
```
1 : 2vqc_A_1_7
2 : 2vqc_A_8_118
```
"""
function _choose_largest(segments)
  j = 0
  len_before = 0
  for i in 1:length(segments)
    len_now = _length_segment(segments[i])
    if len_now > len_before
      j = i
    end
  end
  j
end

"""
Returns an Iterator of the residues on the listResidue
```
<listResidue>
  <residue>
  ...
  </residue>
  ...
</listResidue>
```
"""
_get_residues(segment) = child_elements(select_element(get_elements_by_tagname(segment, "listResidue"), "listResidue"))

"""
Returns `true` if the residue was annotated as *Not_Observed*
"""
function _is_missing(residue)
  details = get_elements_by_tagname(residue, "residueDetail")
  for det in details
    # XML: <residueDetail dbSource="PDBe" property="Annotation">Not_Observed</residueDetail>
    if attribute(det, "property") == "Annotation" && content(det) == "Not_Observed"
      return(true)
    end
  end
  false
end

