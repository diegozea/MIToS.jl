module MSA

  using DataStructures

  export Residue, GAP, @res_str,

  IndexedVector, selectindex, selectvalue, swap!,

  Annotations, filtersequences!, filtercolumns!, empty,
  getannotfile,  getannotcolumn,  getannotsequence,  getannotresidue,
  setannotfile!, setannotcolumn!, setannotsequence!, setannotresidue!,
  annotate_modification!, delete_annotated_modifications!, printmodifications,

  MultipleSequenceAlignment, AnnotatedMultipleSequenceAlignment, AbstractMultipleSequenceAlignment,
  AlignedSequence, AnnotatedAlignedSequence, AbstractAlignedSequence,
  getresidues, getsequence, getresiduesequences,
  nsequences, ncolumns, gappercentage, residuepercentage, coverage,
  columngappercentage, setreference!, gapstrip!, adjustreference!, asciisequence,
  gapstrip, adjustreference, filtersequences, filtercolumns,
  getcolumnmapping, getsequencemapping,

  Raw,

  Stockholm, downloadpfam,

  FASTA,

  shuffle_columnwise!, shuffle_sequencewise!, shuffle_residues_sequencewise!, shuffle_residues_columnwise!

  include("Residues.jl")
  include("IndexedVectors.jl")
  include("Annotations.jl")
  include("MultipleSequenceAlignment.jl")
  include("GeneralParserMethods.jl")
  include("Raw.jl")
  include("Pfam.jl")
  include("FASTA.jl")
  include("Shuffle.jl")

end
