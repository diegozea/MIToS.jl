struct SIFTSXML <: FileFormat end

# Download SIFTS
# ==============

"""
Download the gzipped SIFTS xml  for the `pdbcode`.
The extension of the downloaded file is `.xml.gz` by default.
The `filename` can be changed, but the `.xml.gz` at the end is mandatory.
"""
function downloadsifts(pdbcode::String; filename::String="$(lowercase(pdbcode)).xml.gz")
    @assert endswith(filename, ".xml.gz") "filename must end with .xml.gz"
    if check_pdbcode(pdbcode)
        download_file(string("ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/split_xml/",
            lowercase(pdbcode[2:3]), "/", lowercase(pdbcode), ".xml.gz"), filename)
    else
        throw(ErrorException("$pdbcode is not a correct PDB"))
    end
    filename
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
function _get_entities(sifts)
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
function _get_residues(segment)
    child_elements(select_element(get_elements_by_tagname(segment, "listResidue"), "listResidue"))
end

"""
Returns `true` if the residue was annotated as *Not_Observed*.
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

function _get_details(residue)::Tuple{Bool,String,String}
    details = get_elements_by_tagname(residue, "residueDetail")
    missing_residue = false
    sscode  = " "
    ssname  = " "
    for det in details
        detail_property = attribute(det, "property")
        # XML: <residueDetail dbSource="PDBe" property="Annotation">Not_Observed</residueDetail>
        if  detail_property == "Annotation" && content(det) == "Not_Observed"
            missing_residue = true
            break
        # XML: <residueDetail dbSource="PDBe" property="codeSecondaryStructure"...
        elseif detail_property == "codeSecondaryStructure"
            sscode = content(det)
        # XML: <residueDetail dbSource="PDBe" property="nameSecondaryStructure"...
        elseif detail_property == "nameSecondaryStructure"
            ssname = content(det)
        end
    end
    (missing_residue, sscode, ssname)
end
