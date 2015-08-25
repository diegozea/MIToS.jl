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
    throw(ErrorException(string(pdbcode, " is not a correct PDB")))
  end
end

# Internal Parser Functions
# =========================

"""
Gets the entities of a SIFTS XML. In some cases, each entity is a PDB chain.
WARNING: Sometimes there are more chains than entities!
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

"""
Gets an array of the segments, the continuous region of an entity.
Chimeras and expression tags generates more than one segment for example.
"""
_get_segments(entity) = get_elements_by_tagname(entity, "segment")

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

