using LightXML

function siftsmapping(sifts_xml::ASCIIString, chain::Char, dbid::ASCIIString; db::ASCIIString = "Pfam", segment_num::Int = 1, check_observed::Bool = true)
	mapping = Dict{Int,Int}()
	sifts = parse_file(sifts_xml);
	siftsroot = root(sifts);
	entities = get_elements_by_tagname(siftsroot, "entity");
	for entity in entities
		if attribute(entity, "entityId") == string(chain)
			segments = get_elements_by_tagname(entity, "segment")
			segment_num == 1 && length(segments) != 1 && warn("There is more than one segment, using segment 1: $(attribute(segments[1], "segId"))")
			residues = get_elements_by_tagname(segments[segment_num], "listResidue")
			for residue in child_elements(residues[1])
				use = true
				if check_observed
					details = get_elements_by_tagname(residue, "residueDetail")
					for det in details
						if attribute(det, "property") == "Annotation" && content(det) == "Not_Observed"
							use = false
						end
					end
				end
				if use
					crossref = get_elements_by_tagname(residue, "crossRefDb")
					key = Int( attribute(residue, "dbResNum") ) ## <residue dbSource="PDBe" ...
					for ref in crossref
						if attribute(ref, "dbSource") == db && attribute(ref, "dbAccessionId") == dbid
							mapping[key] = Int( attribute(ref, "dbResNum") )
              break
						end
					end
				end
			end
		end
	end
	sizehint!(mapping, length(mapping))
end

# Download SIFTS
# ==============

function downloadsifts(pdbcode::ASCIIString; filename::ASCIIString="$(lowercase(pdbcode)).xml.gz")
  if length(pdbcode)== 4
    download(string("ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/split_xml/", lowercase(pdbcode[2:3]), "/", lowercase(pdbcode), ".xml.gz"), filename)
  else
    throw(string(pdbcode, " is not a correct PDB"))
  end
end
