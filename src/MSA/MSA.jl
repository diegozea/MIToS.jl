module MSA

  export Residue, GAP, @res_str,

  IndexedVector, selectindex, selectvalue, swap!,

  Annotations, filtersequences!, filtercolumns!, empty,
  getannotfile,  getannotcolumn,  getannotsequence,  getannotresidue,
  setannotfile!, setannotcolumn!, setannotsequence!, setannotresidue!,

  MultipleSequenceAlignment, AnnotatedMultipleSequenceAlignment, AbstractMultipleSequenceAlignment,
  AlignedSequence, AnnotatedAlignedSequence, AbstractAlignedSequence,
  getresidues, getsequence, getresiduesequences,
  nsequences, ncolumns, gappercentage, residuepercentage, coverage,
  columngappercentage, setreference!, gapstrip!, adjustreference!, asciisequence,
  gapstrip, adjustreference, filtersequences, filtercolumns,

  Raw,

  Stockholm, downloadpfam,

  FASTA

  include("Residues.jl")
  include("IndexedVectors.jl")
  include("Annotations.jl")
  include("MultipleSequenceAlignment.jl")
  include("GeneralParserMethods.jl")
  include("Raw.jl")
  include("Pfam.jl")
  include("FASTA.jl")

end
