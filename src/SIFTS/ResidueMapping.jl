using LightXML

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
segment_num = "largest" chooses the largest segment.
segment_num = "iterate" iterates over all the segments.
segment_num = "1" select the segment 1.
"""
function siftsmapping(sifts_xml::ASCIIString, chain::Char, dbid::ASCIIString; db::ASCIIString = "Pfam", segment_num::ASCIIString = "iterate", check_observed::Bool = true)
	mapping = Dict{Int,Int}()
	sifts = parse_file(sifts_xml);
	siftsroot = root(sifts);
	entities = get_elements_by_tagname(siftsroot, "entity");
	for entity in entities
		if attribute(entity, "entityId") == string(chain)
			segments = get_elements_by_tagname(entity, "segment")
      # Choose the segments
      if segment_num == "largest"
        n = _choose_largest(segments)
        indexes = n:n
      elseif segment_num == "iterate"
        indexes = 1:length(segments)
      else
        n = parse(Int, segment_num)
        indexes = n:n
      end
      # Iterates over the selected segments
      for index in indexes

        residues = get_elements_by_tagname(segments[index], "listResidue")
			  for residue in child_elements(residues[1])
				  use = true
				  if check_observed
            # XML: <residueDetail dbSource="PDBe" property="Annotation">Not_Observed</residueDetail>
					  details = get_elements_by_tagname(residue, "residueDetail")
					  for det in details
					  	if attribute(det, "property") == "Annotation" && content(det) == "Not_Observed"
					  		use = false
					  	end
					  end
				  end
				  if use
					  crossref = get_elements_by_tagname(residue, "crossRefDb")
					  key = parse(Int, attribute(residue, "dbResNum")) # XML: <residue dbSource="PDBe" ...
					  for ref in crossref
						  if attribute(ref, "dbSource") == db && attribute(ref, "dbAccessionId") == dbid
							  mapping[key] = parse(Int, attribute(ref, "dbResNum"))
                break
						  end
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
