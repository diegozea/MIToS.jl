module MSA

  export Residue, residue, GAP,

  IndexedVector, indexedvector, selectindex, selectvalue, swap!,

  Annotations, filtersequences!, filtercolumns!,

  MultipleSequenceAlignment, AlignedSequence, getresidues, getsequence,
  nsequences, ncolumns, gappercentage, residuepercentage, coverage, 
  columngappercentage, setreference!, gapstrip!, adjustreference!, asciisequence,

  readpfam, writepfam, printpfam

  include("Residues.jl")
  include("IndexedVectors.jl")
  include("Annotations.jl")
  include("MultipleSequenceAlignment.jl")
  include("Pfam.jl")

end
