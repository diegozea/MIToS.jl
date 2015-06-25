module MSA

  export Residue, residue, GAP,

  IndexedVector, indexedvector, selectindex, selectvalue, swap!,

  Annotations, filtersequences!, filtercolumns!,

  MultipleSequenceAlignment, AlignedSequence, getresidues, getsequence, getrawsequences,
  nsequences, ncolumns, gappercentage, residuepercentage, coverage,
  columngappercentage, setreference!, gapstrip!, adjustreference!, asciisequence,

  readpfam, writepfam, printpfam, downloadpfam,

  readfasta, writefasta, printfasta

  include("Residues.jl")
  include("IndexedVectors.jl")
  include("Annotations.jl")
  include("MultipleSequenceAlignment.jl")
  include("Pfam.jl")
  include("FASTA.jl")

end
