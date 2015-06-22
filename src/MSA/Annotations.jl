# Annotations
# ===========

type Annotations
  file::Dict{ASCIIString, ASCIIString}
  sequences::Dict{(ASCIIString,ASCIIString),ASCIIString}
  columns::Dict{ASCIIString,ASCIIString}
  residues::Dict{(ASCIIString,ASCIIString),ASCIIString}
end

# Filters
# -------

__filter(str::ASCIIString, mask) = ascii( str.data[mask] )

function filtersequences!(data::Annotations, ids::IndexedVector, mask)
  if length(data.sequences) > 0 || length(data.residues) > 0
    del = ids.values[ !mask ]
    for key in keys(data.residues)
      if key[1] in del
        delete!(data.residues, key)
      end
    end
    for key in keys(data.sequences)
      if key[1] in del
        delete!(data.sequences, key)
      end
    end
    data.sequences = sizehint(data.sequences, length(data.sequences))
    data.residues = sizehint(data.residues, length(data.residues))
  end
  data
end

function filtercolumns!(data::Annotations, mask)
  if length(data.columns) > 0 || length(data.residues) > 0
    for (key,value) in data.residues
      data.residues[key] = __filter(value, mask)
    end
    for (key,value) in data.columns
      data.columns[key] = __filter(value, mask)
    end
  end
  data
end

# Copy and deepcopy
# -----------------

deepcopy(ann::Annotations) = Annotations( deepcopy( ann.file ), deepcopy( ann.sequences ), deepcopy( ann.columns ), deepcopy( ann.residues ))

copy(ann::Annotations) = Annotations( copy( ann.file ), copy( ann.sequences ), copy( ann.columns ), copy( ann.residues ))

# Show & Print
# ------------

import Base: print, show

print(io::IO, ann::Annotations) = dump(io, ann)
show(io::IO, ann::Annotations) = dump(io, ann)