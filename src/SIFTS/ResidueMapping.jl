function getsiftsmapping(filename::ASCIIString, chain::Char; from = "PDB", to = "Pfam", segment_num = 1, check_observed = true)
	mapping = Dict{Int,Int}()
	sifts = parse_file(filename);
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
					key = -9999
					value = -9999
					for ref in crossref
						if attribute(ref, "dbSource") == from
							key = int( attribute(ref, "dbResNum") )
						end
						if attribute(ref, "dbSource") == to
							value = int( attribute(ref, "dbResNum") )
						end
						if key != -9999 && value != -9999
							mapping[key] = value
							break
						end
					end
				end
			end
		end
	end
	sizehint(mapping, length(mapping))
end