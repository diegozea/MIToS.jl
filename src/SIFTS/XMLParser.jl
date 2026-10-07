struct SIFTSXML <: FileFormat end

# Download SIFTS
# ==============

function _sifts_url(pdbcode::AbstractString, source::AbstractString)
    @assert source == "ftp" || source == "https" "source must be ftp or https"
    code = Utils._legacy_pdbcode(pdbcode)
    code === nothing && throw(
        ArgumentError(
            "SIFTS XML downloads currently require a legacy PDB ID; $pdbcode has no legacy alias. " *
            "Read a local SIFTSXML file if one is available.",
        ),
    )
    if source == "ftp"
        string(
            "https://ftp.ebi.ac.uk/pub/databases/msd/sifts/split_xml/",
            code[2:3],
            "/",
            code,
            ".xml.gz",
        )
    else
        string("https://www.ebi.ac.uk/pdbe/files/sifts/", code, ".xml.gz")
    end
end

"""
    downloadsifts(pdbcode::AbstractString; filename::AbstractString, source::AbstractString="ftp")

Download the gzipped SIFTS XML file for the provided `pdbcode`.
The downloaded file will have the default extension `.xml.gz`.
While you can change the `filename`, it must include the `.xml.gz` ending.
The `source` keyword argument is set to `"ftp"` by default, downloading from the HTTPS
mirror at `https://ftp.ebi.ac.uk/pub/databases/msd/sifts/split_xml/`.
Alternatively, you can choose `"https"` as the `source` to download directly from the
EBI PDBe server at https://www.ebi.ac.uk/pdbe/files/sifts/.
Legacy IDs and their extended aliases are accepted in either case. The remote SIFTS
XML services still use legacy IDs: `pdb_00001abc` downloads `1abc.xml.gz`, while keeping
`pdb_00001abc.xml.gz` as the default local filename. An extended-only ID raises
`ArgumentError` because no SIFTS XML endpoint for these IDs is currently documented.
"""
function downloadsifts(
    pdbcode::AbstractString;
    filename::AbstractString = "$(lowercase(pdbcode)).xml.gz",
    source::AbstractString = "ftp",
)
    @assert endswith(filename, ".xml.gz") "filename must end with .xml.gz"
    url = _sifts_url(pdbcode, source)
    download_file(url, filename)
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
