## hcat (horizontal concatenation)

function _hcat_seqnames_kernel!(concatenated_names, msa, delim)
	for (i, seq) in enumerate(sequencename_iterator(msa))
		if seq == concatenated_names[i]
			continue
		end
		concatenated_names[i] *= delim * seq
	end
	concatenated_names
end

function _h_concatenated_seq_names(msas...; delim::String="_&_")
	concatenated_names = sequencenames(msas[1])
	for msa in msas[2:end]
		_hcat_seqnames_kernel!(concatenated_names, msa, delim)
	end
	concatenated_names
end

function _hcat_colnames_kernel!(colnames, columns, msa_number::Int)::Int
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
	msa_number
end

function _h_concatenated_col_names(msas...)
	colnames = String[]
	msa_number = 0
	for msa in msas
		columns = columnname_iterator(msa)
		msa_number = _hcat_colnames_kernel!(colnames, columns, msa_number)
	end
	colnames
end

_get_seq_lengths(msas...) = Int[ncolumns(msa) for msa in msas]

function _h_concatenate_annotfile(data::Annotations...)
	N = length(data)
	annotfile = copy(getannotfile(data[1]))
	for i in 2:N
		ann = data[i]::Annotations
		for (k, v) in getannotfile(ann)
			if haskey(annotfile, k)
				if endswith(k, "ColMap") # to also use ColMap annotations from vcat
					annotfile[k] = string(annotfile[k], ",", v)
				else
					annotfile[k] = string(annotfile[k], "_&_", v)
				end
			else
				if endswith(k, "ColMap")
					# that means that the ColMap annotation probably comes from vcat
					# in one of the source MSAs and there is no match between keys.
					@warn "There was no possible to match the ColMap annotations."
				end
				push!(annotfile, k => v)
			end
		end
	end
	annotfile
end

"""
It returns a Dict mapping the MSA number and sequence name to the horizontally 
concatenated sequence name.
"""
function _get_seqname_mapping_hcat(concatenated_seqnames, msas...)
	mapping = Dict{Tuple{Int, String}, String}()
	nmsa = length(msas)
	nseq = length(concatenated_seqnames)
	for j in 1:nmsa
		seq_names = sequencenames(msas[j])
		@assert nseq == length(seq_names)
		for i in 1:nseq
			mapping[(j, seq_names[i])] = concatenated_seqnames[i]
		end
	end
	mapping
end

"""
At this moment, the _concatenate_annotsequence function does no check for SeqMap 
annotations. Therefore, if an MSA does not have a SeqMap annotation, the concatenated
SeqMap annotation is corrupted. To avoid this, we will delete the SeqMap annotation
whose length is not equal to the length of the concatenated MSA.
"""
function _clean_sequence_mapping!(msa::AnnotatedAlignedObject)
	N = ncolumns(msa)
	to_delete = Tuple{String,String}[]
	annotsequence = getannotsequence(msa)
	for (key, value) in annotsequence
		if key[2] == "SeqMap"
			n_delim = count(==(','), value)
			if n_delim != N - 1
				push!(to_delete, key)
			end
		end
	end
	for key in to_delete
		delete!(annotsequence, key)
	end
	msa
end

function _concatenate_annotsequence(seqname_mapping, data::Annotations...)
	annotsequence = Dict{Tuple{String,String},String}()
	for (i, ann::Annotations) in enumerate(data)
		for ((seqname, annot_name), value) in getannotsequence(ann)
			concatenated_seqname = get(seqname_mapping, (i, seqname), seqname)
			new_key = (concatenated_seqname, annot_name)
			# if we used :vcat, new_key will not be present in the dict as the
			# sequence names are disambiguated first
			if haskey(annotsequence, new_key)
				# so, we execute the following code only if we used :hcat
				if annot_name == "SeqMap"
					sep = ","
				else
					sep = "_&_"
				end
				annotsequence[new_key] = string(annotsequence[new_key], sep, value)
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
			dict[key] = string(dict[key], value)
		else
			previous = sum(seq_lengths[(last[key]+1):(i-1)])
			dict[key] = string(dict[key], repeat(" ", previous), value)
		end
		last[key] = i
	else
		if i == 1
			push!(dict, key => value)
		else
			previous = sum(seq_lengths[1:(i-1)])
			push!(dict, key => string(repeat(" ", previous), value))
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

function _h_concatenate_annotcolumn(seq_lengths::Vector{Int}, 
		data::Annotations...)::Dict{String,String}
	annotcolumn = Dict{String,String}()
	last = Dict{String,Int}()
	for (i, ann::Annotations) in enumerate(data)
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

function _delete_hcat_annotfile!(annot::Annotations)
	if haskey(annot.file, "HCat")
		delete!(annot.file, "HCat")
	end
	annot
end

function _set_hcat_annotfile!(annot::Annotations, colnames)
	_delete_hcat_annotfile!(annot)
	setannotfile!(annot, "HCat", 
		join((replace(col, r"_[0-9]+$" => "") for col in colnames), ','))
	annot
end

function Base.hcat(msa::T...) where T <: AnnotatedAlignedObject
	seqnames = _h_concatenated_seq_names(msa...)
	colnames = _h_concatenated_col_names(msa...)
	concatenated_matrix = hcat(getresidues.(msa)...)
	concatenated_msa = _namedresiduematrix(concatenated_matrix, seqnames, colnames)
	seqname_mapping = _get_seqname_mapping_hcat(seqnames, msa...)
	seq_lengths = _get_seq_lengths(msa...)
	old_annot = annotations.([msa...])
	new_annot = Annotations(
		_h_concatenate_annotfile(old_annot...),
		_concatenate_annotsequence(seqname_mapping, old_annot...),
		_h_concatenate_annotcolumn(seq_lengths, old_annot...),
		_h_concatenate_annotresidue(seq_lengths, seqname_mapping, old_annot...)
	)
	_set_hcat_annotfile!(new_annot, colnames)
	new_msa = T(concatenated_msa, new_annot)
	_clean_sequence_mapping!(new_msa)
end

function Base.hcat(msa::T...) where T <: UnannotatedAlignedObject
	concatenated_matrix = hcat(getresidues.(msa)...)
	seqnames = _h_concatenated_seq_names(msa...)
	colnames = _h_concatenated_col_names(msa...)
	concatenated_msa = _namedresiduematrix(concatenated_matrix, seqnames, colnames)
	T(concatenated_msa)
end


"""
It returns a vector of numbers from `1` to N for each column that indicates the source MSA.
The mapping is annotated in the `"HCat"` file annotation of an
`AnnotatedMultipleSequenceAlignment` or in the column names of an `NamedArray` or
`MultipleSequenceAlignment`.

NOTE: When the MSA results from vertically concatenating MSAs using `vcat`, 
the `"HCat"` annotations from the constituent MSAs are renamed as `"1_HCat"`, `"2_HCat"`, 
etc. In that case, the MSA numbers referenced in the column names are provided. 
To access the original annotations, utilize the `getannotfile` function.
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
	colnames = columnname_iterator(msa)
	if !isempty(colnames)
		if !occursin('_', first(colnames))
			throw(ErrorException("Column names have not been generated by `hcat`."))
		end
    	Int[ parse(Int, replace(col, r"_[0-9]+$" => "")) for col in colnames]
	else
		throw(ErrorException("There are not column names!"))
	end
end

gethcatmapping(msa::MultipleSequenceAlignment) = gethcatmapping(namedmatrix(msa))

## vcat (vertical concatenation)

"""
If returns a vector of sequence names for the vertically concatenated MSA. The prefix 
is the number associated to the source MSA. If the sequence name has already a number
as prefix, the MSA number is increased accordingly.
"""
function _v_concatenated_seq_names(msas...; fill_mapping::Bool=false)
	label_mapping = Dict{String,Int}()
	concatenated_names = String[]
	msa_number = 0
	previous_msa_number = 0
	for msa in msas
		msa_label = ""
		msa_number += 1
		for seqname in sequencename_iterator(msa)
			m = match(r"^([0-9]+)_(.*)$", seqname)
			if m === nothing
				msa_label = ""
				new_seqname = "$(msa_number)_$seqname"
			else
				# if the sequence name has already a number as prefix, we increase the
				# MSA number every time the prefix number changes
				current_msa_label = string(m.captures[1])
				if current_msa_label == msa_label
					new_seqname = "$(msa_number)_$(m.captures[2])"
				else
					# avoid increasing the MSA number two times in a row the first time 
					# we find a sequence name with a number as prefix
					if msa_label != ""
						msa_number += 1
					end
					msa_label = current_msa_label
					fill_mapping && push!(label_mapping, msa_label => msa_number)
					new_seqname = "$(msa_number)_$(m.captures[2])"
				end
			end
			previous_msa_number = msa_number
			push!(concatenated_names, new_seqname)
		end
	end
	concatenated_names, label_mapping
end

"""
It returns a Dict mapping the MSA number and sequence name to the vertically 
concatenated sequence name.
"""
function _get_seqname_mapping_vcat(concatenated_seqnames, msas...)
	mapping = Dict{Tuple{Int, String}, String}()
	sequence_number = 0
	for (i, msa) in enumerate(msas)
		for seq in sequencename_iterator(msa)
			sequence_number += 1
			mapping[(i, seq)] = concatenated_seqnames[sequence_number]
		end
	end
	mapping
end

function _update_annotation_name(annot_name, msa_number, label_mapping)
	m = match(r"^([0-9]+)_(.*)$", annot_name)
	if m !== nothing
		# The annotation name has already a number as prefix, so we use the mapping
		# to determine the corresponding MSA number
		if haskey(label_mapping, m.captures[1]) 
			# we avoid taking the MSA number from the name if it is not in the mapping
			# to avoid taking a prefix that was alredy there but related to vcat
			msa_number = label_mapping[m.captures[1]]
		end
		new_annot_name = "$(msa_number)_$(m.captures[2])"
	else
		new_annot_name = "$(msa_number)_$annot_name"
	end
	msa_number, new_annot_name
end

function _v_concatenate_annotfile(label_mapping::Dict{String,Int}, data::Annotations...)
	annotfile = OrderedDict{String,String}()
	msa_number = 0
	for ann::Annotations in data
		msa_number += 1
		for (name, annotation) in getannotfile(ann)
			msa_number, new_name = _update_annotation_name(name, msa_number, label_mapping)
			push!(annotfile, new_name => annotation)
		end
	end
	annotfile
end

"""
Column annotations are disambiguated by adding a prefix to the annotation name as
we do for the sequence names.
"""
function _v_concatenate_annotcolumn(label_mapping::Dict{String,Int}, data::Annotations...)
	annotcolumn = Dict{String,String}()
	msa_number = 0
	for ann::Annotations in data
		msa_number += 1
		for (name, annotation) in getannotcolumn(ann)
			msa_number, new_name = _update_annotation_name(name, msa_number, label_mapping)
			push!(annotcolumn, new_name => annotation)
		end
	end
	annotcolumn
end

"""
Residue annotations are disambiguated by adding a prefix to the sequence name holding
the annotation as we do for the sequence names.
"""
function _v_concatenate_annotresidue(concatenated_seqnames, data::Annotations...)
	annotresidue = Dict{Tuple{String,String},String}()
	for (i, ann::Annotations) in enumerate(data)
		for ((seqname, annot_name), value) in getannotresidue(ann)
			concatenated_seqname = get(concatenated_seqnames, (i, seqname), seqname)
			new_key = (concatenated_seqname, annot_name)
			push!(annotresidue, new_key => value)
		end
	end
	annotresidue
end

function Base.vcat(msa::T...) where T <: AnnotatedAlignedObject
	seqnames, label_mapping = _v_concatenated_seq_names(msa...; fill_mapping=true)
	colnames = columnname_iterator(msa[1])
	concatenated_matrix = vcat(getresidues.(msa)...)
	concatenated_msa = _namedresiduematrix(concatenated_matrix, seqnames, colnames)
	seqname_mapping = _get_seqname_mapping_vcat(seqnames, msa...)
	old_annot = annotations.([msa...])
	new_annot = Annotations(
		_v_concatenate_annotfile(label_mapping, old_annot...),
		_concatenate_annotsequence(seqname_mapping, old_annot...),
		_v_concatenate_annotcolumn(label_mapping, old_annot...),
		_v_concatenate_annotresidue(seqname_mapping, old_annot...)
	)
	# There is no need for a VCat annotation as the source MSA number is already
	# annotated as a prefix in the sequence names
	T(concatenated_msa, new_annot)
end

function Base.vcat(msa::T...) where T <: UnannotatedAlignedObject
	concatenated_matrix = vcat(getresidues.(msa)...)
	seqnames, _ = _v_concatenated_seq_names(msa...)
	colnames = columnname_iterator(msa[1])
	concatenated_msa = _namedresiduematrix(concatenated_matrix, seqnames, colnames)
	T(concatenated_msa)
end

#
# Functions to join, merge or pair MSAs
#
# Currently, these functions will utilize hcat and vcat functions as much as possible.
# If this approach proves to be too slow, we may consider preallocating the result matrices.

# helper functions to name columns in gap blocks
function _get_last_gap_number(name_iterator)
	numbers = (parse(Int, replace(name, "gap:" => "")) 
		for name in name_iterator if startswith(name, "gap:"))
	# as `maximum(numbers, init=0)` doesn't work in Julia versions < 1.6
	# we use `foldl` instead
	foldl(max, numbers, init=0)
end

function _gapblock_columnnames(msa, ncol)
	last_gap_number = _get_last_gap_number(columnname_iterator(msa))
	["gap:$(i)" for i in (last_gap_number+1):(last_gap_number+ncol)]
end

# Since gaps are used for padding, the following function creates an MSA that contains
# only gaps but the proper annotations and names to be used in hcat and vcat.
function _gap_columns(msa, ncol)
	nseq = nsequences(msa)
	empty_mapping = repeat(",", ncol - 1)
	matrix = fill(GAP, nseq, ncol)
	seqnames = sequencenames(msa)
	colnames = _gapblock_columnnames(msa, ncol)
	named_matrix = _namedresiduematrix(matrix, seqnames, colnames)
	block = AnnotatedMultipleSequenceAlignment(named_matrix)
	for seq in seqnames
		setannotsequence!(block, seq, "SeqMap", empty_mapping)
	end
	setannotfile!(block, "ColMap", empty_mapping)
	block
end

function _disambiguate_sequence_names(msa, seqnames)
	disambiguous_seqnames = deepcopy(seqnames)
	last_gap_number = _get_last_gap_number(sequencename_iterator(msa))
	for (i, seqname) in enumerate(seqnames)
		# Check for redundancy and rename if necessary
		if seqname in sequencename_iterator(msa)
			last_gap_number += 1
			disambiguous_seqnames[i] = "gap:$(last_gap_number)"
		end
	end
	disambiguous_seqnames
end



"""
    _renumber_sequence_gaps(msa::AnnotatedMultipleSequenceAlignment)

Renumbers the gap sequences in a given multiple sequence alignment (MSA) object and 
returns a new MSA object with updated sequence names.

This function specifically targets sequences within the MSA that are named with the 
prefix "gap:". It assigns a new sequential number to each gap sequence, ensuring a unique 
and ordered naming convention (e.g., "gap:1", "gap:2", etc.). 

Therefore, this function is useful for renumbering gap sequences in an MSA where gap 
blocks have been inserted from the end to the beginning. This insertion order can lead to 
non-sequential numbering of the gap sequences. For example, in the case of two gap blocks, 
where one block contains a single sequence and the other contains two, the sequences 
might originally be numbered as "gap:3", "gap:1", "gap:2". This function renumbers them 
sequentially so that they are numbered as "gap:1", "gap:2", "gap:3".
"""
function _renumber_sequence_gaps(msa::AnnotatedMultipleSequenceAlignment)
	seqnames = collect(sequencename_iterator(msa))
	gap_number = 0
	old2new = Dict{String,String}()
	for (i, old_seqname) in enumerate(seqnames)
		if startswith(old_seqname, "gap:")
			gap_number += 1
			new_seqname = "gap:$(gap_number)"
			if old_seqname != new_seqname
				old2new[old_seqname] = new_seqname
				seqnames[i] = new_seqname
			end
		end
	end
	annot = _rename_sequences(annotations(msa), old2new)
	named_matrix = _namedresiduematrix(getresidues(msa), seqnames, columnnames(msa))
	AnnotatedMultipleSequenceAlignment(named_matrix, annot)
end

function _renumber_column_gaps(col_names::Vector{String})
	new_colnames = copy(col_names)
	gap_number = 0
	for (i, old_colname) in enumerate(col_names)
		if startswith(old_colname, "gap:")
			gap_number += 1
			new_colnames[i] = "gap:$(gap_number)"
		end
	end
	new_colnames
end

function _renumber_column_gaps(msa::AnnotatedMultipleSequenceAlignment)
	new_colnames = _renumber_column_gaps(columnnames(msa))
	named_matrix = _namedresiduematrix(getresidues(msa), sequencenames(msa), new_colnames)
	AnnotatedMultipleSequenceAlignment(named_matrix, annotations(msa))
end

function _gap_sequences(msa, seqnames)
	nseq = length(seqnames)
	ncol = ncolumns(msa)
	empty_mapping = repeat(",", ncol - 1)
	matrix = fill(GAP, nseq, ncol)
	colnames = columnname_iterator(msa)
	# Names are disambiguated to prevent the "Inconsistent dictionary sizes" error during 
	# vcat. This error arises when the same sequence name appears in different MSAs,
	# as sequence names are utilised as keys in a dictionary.
	disambiguous_seqnames = _disambiguate_sequence_names(msa, seqnames)
	named_matrix = _namedresiduematrix(matrix, disambiguous_seqnames, colnames)
	block = AnnotatedMultipleSequenceAlignment(named_matrix)
	for seq in disambiguous_seqnames
		setannotsequence!(block, string(seq), "SeqMap", empty_mapping)
	end
	block
end

# This functions are used to insert a gap block in a given position in the MSA
# without altering the annotations too much. The idea is the gap blocks are
# no real MSAs, but things inserted into a previously existing MSA.

# Insert gap sequences.

function _get_position(max_pos, position::Int)
	if position < 1
		throw(ArgumentError("The gap block position must be equal or greater than 1."))
	elseif position > max_pos
		max_pos + 1 # to insert the gap block at the end
	else
		position # in the MSA
	end
	position
end

function _vcat_gap_block(msa_a, msa_b)
	matrix_a = getresidues(msa_a)
	matrix_b = getresidues(msa_b)
	concatenated_matrix = vcat(matrix_a, matrix_b)
	seqnames_a = sequencenames(msa_a)
	seqnames_b = sequencenames(msa_b)
	seqnames = vcat(seqnames_a, seqnames_b)
	colnames = columnnames(msa_a)
	named_matrix = _namedresiduematrix(concatenated_matrix, seqnames, colnames)
	annot = merge(annotations(msa_a), annotations(msa_b))
	AnnotatedMultipleSequenceAlignment(named_matrix, annot)
end

function _insert_gap_sequences(msa, seqnames, position)
	gap_block = _gap_sequences(msa, seqnames)
	nseq = nsequences(msa)
	int_position = _get_position(nseq, position)
	if int_position == 1 # at start
		_vcat_gap_block(gap_block, msa)
	elseif int_position == nseq + 1 # after end
		_vcat_gap_block(msa, gap_block)
	else
		_vcat_gap_block(_vcat_gap_block(msa[1:(int_position-1), :], gap_block), 
			msa[int_position:end, :])
	end
end

# Insert gap columns.

function _get_msa_number(colnames, position)
	fields = split(colnames[position], '_')
	if length(fields) == 1
		0 # the column does not have a MSA number as prefix
	else
		parse(Int, first(fields))
	end
end

# rely on hcat to do the job, then correct the annotations
# to avoid changing the MSA index number.
function _fix_msa_numbers(original_msa, int_position, gap_block_columns, gapped_msa)
	# 1. get the MSA number that will be used for the gap block
	original_msa_column_names = columnnames(original_msa)
	ncol = length(original_msa_column_names)
	gap_block_colnames = ["gap:$i" for i in 1:gap_block_columns]
	# 3. conserve the annotations outside the gap block
	new_colnames = if int_position == 1 # at start
		vcat(gap_block_colnames, original_msa_column_names)
	elseif int_position == ncol + 1 # after end
		vcat(original_msa_column_names, gap_block_colnames)
	else
		vcat(original_msa_column_names[1:(int_position-1)], gap_block_colnames, 
			original_msa_column_names[int_position:end])
	end
	# 4. update the names and annotations
	# NOTE: We call _renumber_column_gaps to avoid the problem of duplicated names where
	# more than one gap block is inserted
	renamed_colnames = _renumber_column_gaps(new_colnames)
	setnames!(namedmatrix(gapped_msa), renamed_colnames, 2)
	prev_file_annotations = annotations(original_msa).file
	new_file_annotations = annotations(gapped_msa).file
	if haskey(prev_file_annotations, "HCat")
		_set_hcat_annotfile!(new_file_annotations, renamed_colnames)
	else
		# do not add the HCat annotation if it was not present in the original MSA
		delete!(new_file_annotations, "HCat")
	end
	for key in keys(new_file_annotations)
		if startswith(key, "MIToS_") && !haskey(prev_file_annotations, key)
			# delete new MIToS annotated modifications
			delete!(new_file_annotations, key)
		end
	end
	if haskey(prev_file_annotations, "NCol")
		# do not change the NCol annotation
		new_file_annotations["NCol"] = prev_file_annotations["NCol"]
	end
	gapped_msa
end

function _insert_gap_columns(input_msa, gap_block_columns::Int, position)
	@assert gap_block_columns ≥ 0 "The number of gap columns must be greater than 0."
	msa = AnnotatedMultipleSequenceAlignment(input_msa)
	gap_block_columns == 0 && return msa
	gap_block = _gap_columns(msa, gap_block_columns)
	ncol = ncolumns(msa)
	int_position = _get_position(ncol, position)
	gapped_msa = if int_position == 1 # at start
		hcat(gap_block, msa)
	elseif int_position == ncol + 1 # after end
		hcat(msa, gap_block)
	else
		hcat(msa[:, 1:(int_position -1)], gap_block, msa[:, int_position:end])
	end
	_fix_msa_numbers(msa, int_position, gap_block_columns, gapped_msa)
end

# join for MSAs

"""
	_add_gaps_in_b(msa_a, msa_b, positions_a, positions_b, axis::Int=1)

This is an helper function for the `:left` and `:right` joins. It aligns `msa_b` with
`msa_a` by adding gaps to `msa_b`. Only the sequences or columns in `msa_b` that match
the positions in `msa_a` are kept. Gaps are inserted in `msa_b` at positions where `msa_a`
has sequences or columns that are not matched in `msa_b`. Therefore, we can think of
`msa_a` as the reference MSA and `msa_b` as the MSA to be aligned with the reference.

The `positions_a` and `positions_b` arguments define corresponding positions in `msa_a` 
and `msa_b`, respectively, for accurate gap insertion. The function returns a modified 
`msa_b` with gaps inserted to align with `msa_a`.

The axis argument (`1` for sequences, `2` for columns) determines whether gaps are 
inserted as sequences or columns.

# Example
```julia
aligned_msa_b = _add_gaps_in_b(msa_a, msa_b, [1, 2, 3], [2, 3, 4], 2) # axis = 2
```
"""
function _add_gaps_in_b(msa_a, msa_b, positions_a, positions_b, axis::Int=1)
	@assert axis == 1 || axis == 2 "The axis must be 1 (sequences) or 2 (columns)."
	# sort the positions to keep the order in msa_a
	order_a = sortperm(positions_a)
	sorted_b = positions_b[order_a]
	# keep only the matching positions in msa_b
	matching_b = if axis == 1
		AnnotatedMultipleSequenceAlignment(msa_b[sorted_b, :])
	else
		AnnotatedMultipleSequenceAlignment(msa_b[:, sorted_b])
	end
	# find the positions were we need to add gaps in msa_b
	N_a = axis == 1 ? nsequences(msa_a) : ncolumns(msa_a)
	names_a = axis == 1 ? sequencenames(msa_a) : columnnames(msa_a)
	positions_a_set = Set{Int}(positions_a)
	last_b = 0
	gap_positions = Int[]
	gap_names = String[]
	for i in 1:N_a
		if i in positions_a_set # we have seen a matched position in msa_b
			last_b += 1
		else # if the msa_a position is not matched, add gaps in msa_b
			push!(gap_positions, last_b)
			push!(gap_names, names_a[i])
		end
	end
	# compress the list of gap positions to later add full gap blocks
	for pos in Iterators.reverse(unique(gap_positions))
		block_names = gap_names[gap_positions .== pos]
		if axis == 1
			matching_b = _renumber_sequence_gaps(
				_insert_gap_sequences(matching_b, block_names, pos+1))
		else
			matching_b = _renumber_column_gaps(
				_insert_gap_columns(matching_b, length(block_names), pos+1))
		end
	end
	matching_b
end

function _find_pairing_positions(index_function::Function, msa_a, msa_b, pairing)
	n = length(pairing)
	positions_a = Vector{Int}(undef, n)
	positions_b = Vector{Int}(undef, n)
	for (i, (a, b)) in enumerate(pairing)
		positions_a[i] = index_function(msa_a, a)
		positions_b[i] = index_function(msa_b, b)
	end
	positions_a, positions_b
end

function _find_pairing_positions(axis::Int, msa_a, msa_b, pairing)
	@assert axis == 1 || axis == 2 "The axis must be 1 (sequences) or 2 (columns)."
	index_function = axis == 1 ? sequence_index : column_index
	_find_pairing_positions(index_function, msa_a, msa_b, pairing)
end

"""
	_find_gaps(positions, n)

Calculate gaps in a sorted sequence of positions given also a maximum value `n` that would
be added at the end as `n+1` if `n` it is not already present.

This function returns a list of tuples, where each tuple represents a gap in the sequence. 
The first element of the tuple indicates the position after which the gap will be inserted, 
and the second element is the position that comes just after the gap. The function accounts 
for gaps at both the beginning and end of the sequence. A gap is identified as a difference 
greater than 1 between consecutive positions. To account for gaps at the end, the end 
position of a gap matching the last position of the sequence is set to `n+1`.
"""
function _find_gaps(positions, n)
	_positions = if positions[end] == n
		((0,), positions)
	else
		((0,), positions, (n + 1,))
	end
	pos_iterator =  Iterators.flatten(_positions)
	[ (x, y) for (x, y) in zip(Iterators.drop(pos_iterator, 1), pos_iterator) if x - y > 1 ]
end

"""
    _insert_sorted_gaps(msa_target, msa_reference, positions_target, positions_reference, 
		block_position::Symbol=:before, axis::Int=1)

Inserts gap blocks into `msa_target` to match the alignment pattern of `msa_reference`. 
This function is utilized for aligning two MSAs based on specific alignment positions. The 
`positions_target` and `positions_reference` parameters dictate the corresponding 
positions in both MSAs for accurate gap insertion.

The function returns the modified `msa_target` with inserted gaps, aligned according 
to `msa_reference`.

The `block_position` keyword argument determines whether the gap blocks are inserted
`:before` (the default) or `:after` the unmatched sequences or columns. The `axis` keyword 
argument determines whether the gap blocks are inserted as sequences (`1`, the default) or 
columns (`2`).

# Example
```julia
gapped_msa_a = _insert_sorted_gaps(msa_a, msa_b, positions_a, positions_b)
```
"""
function _insert_sorted_gaps(msa_target, msa_reference, positions_target, 
		positions_reference; block_position::Symbol=:before, axis::Int=1)
	@assert block_position in [:before, :after] "block_position must be :before or :after."
	@assert axis == 1 || axis == 2 "The axis must be 1 (sequences) or 2 (columns)."
    # Obtain the positions that will be aligned to gaps in the reference
	N_reference = axis == 1 ? nsequences(msa_reference) : ncolumns(msa_reference)
    gaps_reference = _find_gaps(positions_reference, N_reference)
    # We need the sequence names from the reference that will be aligned to gaps
    names_reference = axis == 1 ? sequencenames(msa_reference) : columnnames(msa_reference)
    # Create a dictionary to find the matching position in `msa_target` for adding the gap blocks
    reference2target = Dict{Int,Int}(positions_reference .=> positions_target)
	N_target = axis == 1 ? nsequences(msa_target) : ncolumns(msa_target)
    # Add the gap blocks to `msa_target` as found when looking at `positions_reference`
    blocks_target = Tuple{Int,Vector{String}}[]
    for (stop, start) in gaps_reference
        # Found the matching position in `msa_target`
        start_target = start == 0 ? 1 : reference2target[start] + 1
		stop_target = stop > N_reference ? N_target + 1 : reference2target[stop]
        # This should work fine, even if `start` is 0 and `stop` is n+1
        unmatched_names = names_reference[start+1:stop-1]
		# Determine whether the gap block should be inserted before or after the sequences
		if block_position == :before
        	push!(blocks_target, (start_target, unmatched_names))
		else
			push!(blocks_target, (stop_target, unmatched_names))
		end
	end
    gapped_msa_target = AnnotatedMultipleSequenceAlignment(msa_target)
	if axis == 1
    	for (position, seqnames) in Iterators.reverse(blocks_target)
        	gapped_msa_target = _insert_gap_sequences(gapped_msa_target, seqnames, position)
    	end
		_renumber_sequence_gaps(gapped_msa_target)
	else
		for (position, colnames) in Iterators.reverse(blocks_target)
			gapped_msa_target = _insert_gap_columns(gapped_msa_target, length(colnames), position)
		end
		_renumber_column_gaps(gapped_msa_target)
	end
end

function _reorder_and_extract_unmatched_names(msa, positions, axis::Int)
    N = size(msa, axis)
    unmatched_positions = setdiff(1:N, positions)
    if axis == 1
		reordered_msa = msa[vcat(positions, unmatched_positions), :]
		reordered_names = sequencenames(reordered_msa)
	else
		reordered_msa = msa[:, vcat(positions, unmatched_positions)]
		reordered_names = columnnames(reordered_msa)
	end
    unmatched_names = reordered_names[length(positions)+1:end]

    return reordered_msa, unmatched_names
end

function Base.join(msa_a::AnnotatedMultipleSequenceAlignment, 
	msa_b::AnnotatedMultipleSequenceAlignment, pairing; kind::Symbol=:outer, axis::Int=1)
	# Check that input arguments and throw explicit ArgumentErrors if necessary
	if isempty(pairing)
		throw(ArgumentError("The pairing argument indicating the matching positions is empty."))
	end
	if length(first(pairing)) != 2
		throw(ArgumentError(string("pairing must consist of pairs where the first element ",
		"corresponds to a position in the first MSA and the second to the second MSA. ", 
		"Otherwise, utilize two distinct position lists.")))
	end
	if kind != :inner && kind != :outer && kind != :left && kind != :right
		throw(ArgumentError("The kind of join must be one of :inner, :outer, :left or :right."))
	end
	if axis != 1 && axis != 2
		throw(ArgumentError("The axis must be 1 (sequences) or 2 (columns)."))
	end
	positions_a, positions_b = _find_pairing_positions(axis, msa_a, msa_b, pairing)
	if kind == :inner
		if axis == 1
			return hcat(msa_a[positions_a, :], msa_b[positions_b, :])
		else
			return vcat(msa_a[:, positions_a], msa_b[:, positions_b])
		end
	elseif kind == :outer
		if issorted(positions_a) && issorted(positions_b)
			gapped_a = _insert_sorted_gaps(msa_a, msa_b, positions_a, positions_b, 
				block_position=:after, axis=axis)
			gapped_b = _insert_sorted_gaps(msa_b, msa_a, positions_b, positions_a, 
				axis=axis)
			return axis == 1 ? hcat(gapped_a, gapped_b) : vcat(gapped_a, gapped_b)
		else
			# do as the inner join for the matching positions, and add the unmatched
			# positions at the end, using gap blocks
			reordered_a, unmatched_names_a = _reorder_and_extract_unmatched_names(msa_a, 
				positions_a, axis)
			reordered_b, unmatched_names_b = _reorder_and_extract_unmatched_names(msa_b, 
				positions_b, axis)
			if axis == 1
				a_block = _insert_gap_sequences(reordered_a, unmatched_names_b, 
					nsequences(reordered_a) + 1)
				b_block = _insert_gap_sequences(reordered_b, unmatched_names_a,
					length(positions_b) + 1)
				return hcat(a_block, b_block)
			else
				a_block = _insert_gap_columns(reordered_a, length(unmatched_names_b), 
					ncolumns(reordered_a) + 1)
				b_block = _insert_gap_columns(reordered_b, length(unmatched_names_a),
					length(positions_b) + 1)
				return vcat(a_block, b_block)
			end
		end
	elseif kind == :left
		gapped_b = _add_gaps_in_b(msa_a, msa_b, positions_a, positions_b, axis)
		return axis == 1 ? hcat(msa_a, gapped_b) : vcat(msa_a, gapped_b)
	elseif kind == :right
		gapped_a = _add_gaps_in_b(msa_b, msa_a, positions_b, positions_a, axis)
		return axis == 1 ? hcat(gapped_a, msa_b) : vcat(gapped_a, msa_b)
	else
		throw(ArgumentError("The kind of join must be one of :inner, :outer, :left or :right."))
	end
end

function Base.join(msa_a::AnnotatedMultipleSequenceAlignment, 
		msa_b::AnnotatedMultipleSequenceAlignment, positions_a, positions_b; 
		kind::Symbol=:outer, axis::Int=1)
	if length(positions_a) != length(positions_b)
		throw(ArgumentError("positions_a and positions_b must have the same length."))
	end
	pairing = zip(positions_a, positions_b)
	join(msa_a, msa_b, pairing, kind=kind, axis=axis)
end