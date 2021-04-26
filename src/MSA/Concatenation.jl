function _concatenated_seq_names(msas...; delim::String="_&_")
	concatenated_names = sequencenames(msas[1])
	for msa in msas[2:end]
		for (i, seq) in enumerate(sequencenames(msa))
			if seq == concatenated_names[i]
				continue
			end
			concatenated_names[i] *= delim * seq
		end
	end
	concatenated_names
end

function _concatenated_col_names(msas...)
	colnames = String[]
	for (i, msa) in enumerate(msas)
		for col in columnnames(msa)
			push!(colnames, "$(i)_$col")
		end
	end
	colnames
end

function _concatenate_annotfile(data::Annotations...)
	annotfile = copy(getannotfile(data[1]))
	for ann in data[2:end]
		for (k, v) in getannotfile(ann)
			if haskey(annotfile, k)
				if k == "ColMap"
					annotfile[k] *= ',' * v
				else
					annotfile[k] *= "_&_" * v
				end
			else
				push!(annotfile, k => v)
			end
		end
	end
	annotfile
end

function _get_seqname_mapping(concatenated_seqnames)
	mapping = Dict{Tuple{Int, String}, String}()
	for concatenated_seqname in concatenated_seqnames
		for (i, seqname) in enumerate(split(concatenated_seqname, "_&_"))
			mapping[(i, seqname)] = concatenated_seqname
		end
	end
	mapping
end

function _concatenate_annotsequence(seqname_mapping, data::Annotations...)
	annotsequence = Dict{Tuple{String,String},String}()
	for (i, ann) in enumerate(data)
		for ((seqname, annot_name), value) in getannotsequence(ann)
			if haskey(seqname_mapping, (i, seqname))
				concatenated_seqname = seqname_mapping[(i, seqname)]
				new_key = (concatenated_seqname, annot_name)
				if haskey(annotsequence, new_key)
					if annot_name == "SeqMap"
						annotsequence[new_key] *= ',' * value
					else
						annotsequence[new_key] *= "_&_" * value
					end
				else
					push!(annotsequence, new_key => value)
				end
			end
		end
	end
	annotsequence
end

function _concatenate_annotcolumn(data::Annotations...)
	annotcolumn = copy(getannotcolumn(data[1]))
	for ann in data[2:end]
		for (k, v) in getannotcolumn(ann)
			if haskey(annotcolumn, k)
				annotcolumn[k] *= v
			else
				push!(annotcolumn, k => v)
			end
		end
	end
	annotcolumn
end

function _concatenate_annotresidue(seqname_mapping, data::Annotations...)
	annotresidue = Dict{Tuple{String,String},String}()
	for (i, ann) in enumerate(data)
		for ((seqname, annot_name), value) in getannotresidue(ann)
			if haskey(seqname_mapping, (i, seqname))
				concatenated_seqname = seqname_mapping[(i, seqname)]
				new_key = (concatenated_seqname, annot_name)
				if haskey(annotresidue, new_key)
					annotresidue[new_key] *= value
				else
					push!(annotresidue, new_key => value)
				end
			end
		end
	end
	annotresidue
end

function Base.hcat(msa::T...) where T <: AnnotatedAlignedObject
	concatenated_msa = NamedArray(hcat(getresidues.(msa)...))
	seqnames = _concatenated_seq_names(msa...)
	colnames = _concatenated_col_names(msa...)
	setnames!(concatenated_msa, seqnames, 1)
	setnames!(concatenated_msa, colnames, 2)
	seqname_mapping = _get_seqname_mapping(seqnames)
	old_annot = annotations.([msa...])
	new_annot = Annotations(
		_concatenate_annotfile(old_annot...),
		_concatenate_annotsequence(seqname_mapping, old_annot...),
		_concatenate_annotcolumn(old_annot...),
		_concatenate_annotresidue(seqname_mapping, old_annot...)
	)
	setannotcolumn!(
		new_annot, 
		"HCat", 
		join(replace(col, r"_[0-9]+$" => "") for col in colnames)
	)
	T(concatenated_msa, new_annot)
end

function Base.hcat(msa::T...) where T <: UnannotatedAlignedObject
	concatenated_msa = NamedArray(hcat(getresidues.(msa)...))
	seqnames = _concatenated_seq_names(msa...)
	colnames = _concatenated_col_names(msa...)
	setnames!(concatenated_msa, seqnames, 1)
	setnames!(concatenated_msa, colnames, 2)
	T(concatenated_msa)
end
