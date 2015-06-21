module MSA

  export Residue, residue, GAP,

  IndexedVector, indexedvector, selectindex, selectvalue, swap!,

  Annotations, filtersequences!, filtercolumns!,

  MultipleSequenceAlignment, AlignedSequence, getresidues,
  nsequences, ncolumns, getsequence, gappercentage, coverage,
  columngappercentage, setreference!, gapstrip!, adjustreference!,

  readpfam

  include("Residues.jl")
  include("IndexedVectors.jl")
  include("Annotations.jl")
  include("MultipleSequenceAlignment.jl")
  include("Pfam.jl")

end
