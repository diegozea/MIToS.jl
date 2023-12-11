## hcat (horizontal concatenation)

function _h_concatenated_seq_names(msas...; delim::String="_&_")
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

function _h_concatenated_col_names(msas...)
	colnames = String[]
	msa_number = 0
	for msa in msas
		columns = columnnames(msa)
		first_col = first(columns)
		check_msa_change = '_' in first_col
		previous = ""
		msa_number += 1
		for col in columns
			if check_msa_change
				fields = split(col, '_')
				current = first(fields)
				if current != previous
					if previous != ""
						msa_number += 1
					end
					previous = string(current)
				end
				col = last(fields)
			end
			push!(colnames, "$(msa_number)_$col")
		end
	end
	colnames
end

_get_seq_lengths(msas...) = [ncolumns(msa) for msa in msas]

function _get_annot_types(fun, index, data::Annotations...)
	union(Set(k[index] for k in keys(fun(annot))) for annot in data)
end

function check_mode(mode::Symbol)
	if mode != :hcat && mode != :vcat
		throw(ArgumentError("The mode must be :hcat or :vcat"))
	end
	mode
end

function _concatenate_annotfile(data::Annotations...; mode::Symbol=:hcat)
	check_mode(mode)
	annotfile = copy(getannotfile(data[1]))
	for ann in data[2:end]
		for (k, v) in getannotfile(ann)
			if haskey(annotfile, k)
				if k == "ColMap"
					if mode == :hcat
						annotfile[k] *= ',' * v
					end
					# if the mode is :vcat, annotfile[k] is not modified, so that the first
					# ColMap is kept
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

"""
It returns a Dict mapping the MSA number and sequence name to the concatenated sequence name.
"""
function _get_seqname_mapping(concatenated_seqnames, msas...)
	mapping = Dict{Tuple{Int, String}, String}()
	seq_names = hcat([sequencenames(msa) for msa in msas]...)
	nseq, nmsa = size(seq_names)
	@assert nseq == length(concatenated_seqnames)
	for i in 1:nseq
		for j in 1:nmsa
			mapping[(j, seq_names[i, j])] = concatenated_seqnames[i]
		end
	end
	mapping
end

function _h_concatenate_annotsequence(seqname_mapping, data::Annotations...)
	annotsequence = Dict{Tuple{String,String},String}()
	for (i, ann) in enumerate(data)
		for ((seqname, annot_name), value) in getannotsequence(ann)
			concatenated_seqname = get(seqname_mapping, (i, seqname), seqname)
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
	annotsequence
end

function _fill_and_update!(dict, last, key, i, value, seq_lengths)
	if haskey(dict, key)
		if last[key] == i - 1
			dict[key] *= value
		else
			previous = sum(seq_lengths[(last[key]+1):(i-1)])
			dict[key] *= " "^previous * value
		end
		last[key] = i
	else
		if i == 1
			push!(dict, key => value)
		else
			previous = sum(seq_lengths[1:(i-1)])
			push!(dict, key => " "^previous * value)
		end
		push!(last, key => i)
	end
end

function _fill_end!(dict, seq_lengths, entity)
	total_length = sum(seq_lengths)
	for (key, value) in dict
		current_length = length(value)
		if current_length < total_length
			dict[key] *= " "^(total_length - current_length)
		elseif current_length > total_length
			throw(ErrorException(
				"There are $current_length $entity annotated instead of $total_length."))
		end
	end
	dict
end

function _h_concatenate_annotcolumn(seq_lengths, data::Annotations...)
	annotcolumn = Dict{String,String}()
	last = Dict{String,Int}()
	for (i, ann) in enumerate(data)
		for (annot_name, value) in getannotcolumn(ann)
			_fill_and_update!(annotcolumn, last, annot_name, i, value, seq_lengths)
		end
	end
	_fill_end!(annotcolumn, seq_lengths, "columns")
end

function _h_concatenate_annotresidue(seq_lengths, seqname_mapping, data::Annotations...)
	annotresidue = Dict{Tuple{String,String},String}()
	last = Dict{Tuple{String,String},Int}()
	for (i, ann) in enumerate(data)
		for ((seqname, annot_name), value) in getannotresidue(ann)
			concatenated_seqname = get(seqname_mapping, (i, seqname), seqname)
			new_key = (concatenated_seqname, annot_name)
			_fill_and_update!(annotresidue, last, new_key, i, value, seq_lengths)
		end
	end
	_fill_end!(annotresidue, seq_lengths, "residues")
end

function Base.hcat(msa::T...) where T <: AnnotatedAlignedObject
	concatenated_msa = NamedArray(hcat(getresidues.(msa)...))
	seqnames = _h_concatenated_seq_names(msa...)
	colnames = _h_concatenated_col_names(msa...)
	setnames!(concatenated_msa, seqnames, 1)
	setnames!(concatenated_msa, colnames, 2)
	seqname_mapping = _get_seqname_mapping(seqnames, msa...)
	seq_lengths = _get_seq_lengths(msa...)
	old_annot = annotations.([msa...])
	new_annot = Annotations(
		_concatenate_annotfile(old_annot...; mode=:hcat),
		_h_concatenate_annotsequence(seqname_mapping, old_annot...),
		_h_concatenate_annotcolumn(seq_lengths, old_annot...),
		_h_concatenate_annotresidue(seq_lengths, seqname_mapping, old_annot...)
	)
	if haskey(new_annot.file, "HCat")
		delete!(new_annot.file, "HCat")
	end
	setannotfile!(
		new_annot, 
		"HCat", 
		join((replace(col, r"_[0-9]+$" => "") for col in colnames), ',')
	)
	T(concatenated_msa, new_annot)
end

function Base.hcat(msa::T...) where T <: UnannotatedAlignedObject
	concatenated_msa = NamedArray(hcat(getresidues.(msa)...))
	seqnames = _h_concatenated_seq_names(msa...)
	colnames = _h_concatenated_col_names(msa...)
	setnames!(concatenated_msa, seqnames, 1)
	setnames!(concatenated_msa, colnames, 2)
	T(concatenated_msa)
end


"""
It returns a vector of numbers from `1` to N for each column that indicates the source MSA.
The mapping is annotated in the `"HCat"` file annotation of an
`AnnotatedMultipleSequenceAlignment` or in the column names of an `NamedArray` or
`MultipleSequenceAlignment`.
"""
function gethcatmapping(msa::AnnotatedMultipleSequenceAlignment)
    annot = getannotfile(msa)
    if haskey(annot, "HCat")
        return _str2int_mapping(annot["HCat"])
    else
        return gethcatmapping(namedmatrix(msa))
    end
end

function gethcatmapping(msa::NamedResidueMatrix{AT}) where AT <: AbstractMatrix
	colnames = names(msa, 2)
	if !isempty(colnames)
		if !occursin('_', colnames[1])
			throw(ErrorException(
				"The column names have not generated by `hcat` on an `AnnotatedMultipleSequenceAlignment` or `MultipleSequenceAlignment`."
				))
		end
    	Int[ parse(Int, replace(col, r"_[0-9]+$" => "")) for col in colnames]
	else
		throw(ErrorException("There are not column names!"))
	end
end

gethcatmapping(msa::MultipleSequenceAlignment) = gethcatmapping(namedmatrix(msa))

## vcat (vertical concatenation)

"""
If returns a vector of sequence names for the vertically concatenated MSA. The names are
kept the same if no sequence name is repeated. If there are repeated sequence names, a 
disambiguation prefix is added to the sequence names. The prefix is the number associated
to the source MSA.
"""
function _v_concatenated_seq_names(msas...)
	seqnames = sequencenames(msas[1])
	msa_number = 1
	for msa in msas[2:end]
		msa_number += 1
		for seq in sequencenames(msa)
			if seq in seqnames
				seq = "$(msa_number)_$seq"
			end
			push!(seqnames, seq)
		end
	end
	seqnames
end



## join



