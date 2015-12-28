module MSA

  using DataStructures
  using IndexedArrays
  using MIToS.Utils

  "`swap!(ia::IndexedArray, to::Int, from::Int)` interchange/swap the values on the indices `to` and `from` in the `IndexedArray`"
  function swap!(ia::IndexedArray, to::Int, from::Int)
    previous_id  = ia[to]
    future_id    = ia[from]

    ia.items[to]   = future_id
    ia.items[from] = previous_id

    ia.lookup[previous_id] = from
    ia.lookup[future_id]   = to

    ia
  end

  import Base: parse, print, write

  export Residue, GAP, @res_str,
  swap!,

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

  Raw, Stockholm, FASTA,

  shuffle_columnwise!, shuffle_sequencewise!, shuffle_residues_sequencewise!, shuffle_residues_columnwise!

  include("Residues.jl")
  include("Annotations.jl")
  include("MultipleSequenceAlignment.jl")
  include("GeneralParserMethods.jl")
  include("Raw.jl")
  include("Stockholm.jl")
  include("FASTA.jl")
  include("Shuffle.jl")

end
