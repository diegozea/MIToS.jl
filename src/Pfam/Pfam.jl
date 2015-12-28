module Pfam

  using MIToS.Utils
  using MIToS.MSA

  export Stockholm, downloadpfam, getseq2pdb

  include("download.jl")
  include("pdb.jl")

end
