struct SIFTSXML <: FileFormat end

# Download SIFTS
# ==============

"""
    downloadsifts(pdbcode::String; filename::String, source::String="https")

Download the gzipped SIFTS XML file for the provided `pdbcode`.
The downloaded file will have the default extension `.xml.gz`.
While you can change the `filename`, it must include the `.xml.gz` ending.
The `source` keyword argument is set to `"https"` by default.
Alternatively, you can choose `"ftp"` as the `source`, which will retrieve the file from
the EBI FTP server at ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/.
However, please note that using `"https"` is highly recommended.
This option will download the file from the
EBI PDBe server at https://www.ebi.ac.uk/pdbe/files/sifts/.
"""
function downloadsifts(
    pdbcode::String;
    filename::String = "$(lowercase(pdbcode)).xml.gz",
    source::String = "https",
)
    @assert endswith(filename, ".xml.gz") "filename must end with .xml.gz"
    @assert source == "ftp" || source == "https" "source must be ftp or https"
    if check_pdbcode(pdbcode)
        url = if source == "ftp"
            string(
                "ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/split_xml/",
                lowercase(pdbcode[2:3]),
                "/",
                lowercase(pdbcode),
                ".xml.gz",
            )
        else
            string("https://www.ebi.ac.uk/pdbe/files/sifts/", lowercase(pdbcode), ".xml.gz")
        end
        download_file(url, filename)
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
    siftsroot = LightXML.root(sifts)
    LightXML.get_elements_by_tagname(siftsroot, "entity")
end

"""
Gets an array of the segments, the continuous region of an entity.
Chimeras and expression tags generates more than one segment for example.
"""
_get_segments(entity) = LightXML.get_elements_by_tagname(entity, "segment")

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
    LightXML.child_elements(
        select_element(
            LightXML.get_elements_by_tagname(segment, "listResidue"),
            "listResidue",
        ),
    )
end

"""
Returns `true` if the residue was annotated as *Not_Observed*.
"""
function _is_missing(residue)
    details = LightXML.get_elements_by_tagname(residue, "residueDetail")
    for det in details
        # XML: <residueDetail dbSource="PDBe" property="Annotation">Not_Observed</residueDetail>
        if LightXML.attribute(det, "property") == "Annotation" &&
           LightXML.content(det) == "Not_Observed"
            return (true)
        end
    end
    false
end

function _get_details(residue)::Tuple{Bool,String,String}
    details = LightXML.get_elements_by_tagname(residue, "residueDetail")
    missing_residue = false
    sscode = " "
    ssname = " "
    for det in details
        detail_property = LightXML.attribute(det, "property")
        # XML: <residueDetail dbSource="PDBe" property="Annotation">Not_Observed</residueDetail>
        if detail_property == "Annotation" && LightXML.content(det) == "Not_Observed"
            missing_residue = true
            break
            # XML: <residueDetail dbSource="PDBe" property="codeSecondaryStructure"...
        elseif detail_property == "codeSecondaryStructure"
            sscode = LightXML.content(det)
            # XML: <residueDetail dbSource="PDBe" property="nameSecondaryStructure"...
        elseif detail_property == "nameSecondaryStructure"
            ssname = LightXML.content(det)
        end
    end
    (missing_residue, sscode, ssname)
end
