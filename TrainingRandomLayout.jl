#!/usr/bin/env julia

import PairwiseListMatrices
@everywhere using PairwiseListMatrices
import DataFrames
@everywhere using DataFrames
import MIToS.Pfam
@everywhere using MIToS.Pfam
import MIToS.MSA
@everywhere using MIToS.MSA
import MIToS.Information
@everywhere using MIToS.Information
import MIToS.Clustering
@everywhere using MIToS.Clustering
import MIToS.SIFTS
@everywhere using MIToS.SIFTS
import MIToS.PDB
@everywhere using MIToS.PDB
import Distributions
@everywhere using Distributions


@everywhere begin
  const data = readtable("Training.csv")

  const pfamdir  = "/data/diego/Pfam"
  const msafiles = filter!(r".*\.aligned\.gz$", readdir(pfamdir))
  const msafile  = [ ascii(split(f, '.')[1]) => joinpath(pfamdir,f) for f in msafiles ];

  const pdbdir   = "/data/diego/PDB"
  const pdbfiles = filter!(r"\S+\.xml\.gz$", readdir(pdbdir))
  const pdbfile  = [ ascii(split(f, '.')[1]) => joinpath(pdbdir,f) for f in pdbfiles ];

  const siftsdir   = "/home/diego/SIFTS"
  const siftsfiles = filter!(r"\S+\.sifts\.xml\.gz$", readdir(siftsdir))
  const siftsfile  = [ uppercase(ascii(split(f, '.')[1])) => joinpath(siftsdir,f) for f in siftsfiles ];

  function _get_msa(pfam_acc, seqid)
    msa = read(msafile[pfam_acc], Stockholm) # generatemapping=true, useidcoordinates=true
    setreference!(msa, seqid)
  end

  function _get_res(pdbcode, chain, pdbdir)
    if haskey(pdbfile, pdbcode)
      pdb = residuesdict(read(pdbfile[pdbcode], PDBML,
                              chain=ascii(chain), model="1", group="ATOM", onlyheavy=true),
                         "1", ascii(chain), "ATOM", "*")
    else
      outfile = joinpath(pdbdir,string(pdbcode,".xml.gz"))
      downloadpdb(pdbcode, outfile=outfile)
      pdb = residuesdict(read(outfile, PDBML,
                              chain=ascii(chain), model="1", group="ATOM", onlyheavy=true),
                         "1", ascii(chain), "ATOM", "*")
    end
    pdb
  end

  function _get_map(msa, seqid, pdbcode, chain, pfam, siftsdir)
    if haskey(siftsfile, pdbcode)
      map = msacolumn2pdbresidue(msa, ascii(seqid), uppercase(pdbcode), ascii(chain), ascii(pfam), ascii(siftsfile[pdbcode]))
    else
      outfile = ascii(string(joinpath(siftsdir,lowercase(pdbcode)),".sifts.xml.gz"))
      downloadsifts(ascii(pdbcode), filename=outfile)
      map = msacolumn2pdbresidue(msa, ascii(seqid), uppercase(pdbcode), ascii(chain), ascii(pfam), outfile)
    end
    map
  end

  function setup!(pfam, seqid, pdb, chain, pdbdir, siftsdir)
    msa = _get_msa(pfam, seqid)
    res = _get_res(pdb, chain, pdbdir)
    map = _get_map(msa, seqid, pdb, chain, pfam, siftsdir)
    pdbmask = hasresidues(msa, map)
    coverage = mean(pdbmask)
    filtercolumns!(msa, pdbmask)
    gapstrip!(msa)
    map = _get_map(msa, seqid, pdb, chain, pfam, siftsdir)
    contactmap = msacontacts(msa, res, map)
    truecontacts, notcontacts = getcontactmasks(contactmap)
    msa, truecontacts, notcontacts, coverage
  end

  function _blmi(msa, beta, clustering)
    aln = copy(msa)
    numbercl = getnclusters(clustering)
    mi = estimateincolumns(aln, Float64, ResidueProbability{Float64, 2, false}, numbercl, beta,
                           MutualInformation{Float64}(), AdditiveSmoothing(0.0), clustering, false, 0.0)
    APC!(mi)
    rand_mi = Array(PairwiseListMatrix{Float64,false}, 20)
    for ns in 1:20
      shuffle_residues_sequencewise!(aln)
      rand_mi[ns] = APC!(estimateincolumns(aln, Float64, ResidueProbability{Float64, 2, false}, numbercl, beta,
                                           MutualInformation{Float64}(), AdditiveSmoothing(0.0), clustering, false, 0.0))
    end
    PairwiseListMatrices.zscore(rand_mi, mi), mi
  end

  function _buslje09(msa, lambda, clustering)
    aln = copy(msa)
    mi = estimateincolumns(aln, ResidueCount{Float64, 2, false}, MutualInformation{Float64}(),
                           AdditiveSmoothing{Float64}(lambda), clustering, false)
    APC!(mi)
    rand_mi = Array(PairwiseListMatrix{Float64,false}, 20)
    for ns in 1:20
      shuffle_residues_sequencewise!(aln)
      rand_mi[ns] = APC!(estimateincolumns(aln, ResidueCount{Float64, 2, false}, MutualInformation{Float64}(),
                                           AdditiveSmoothing{Float64}(lambda), clustering, false))
    end
    PairwiseListMatrices.zscore(rand_mi, mi), mi
  end

  function optim_row(row)

    uniprot = ascii(row[:uniprot])
    chain = ascii(row[:chain])
    pdb = ascii(row[:pdb])
    pfam = ascii(row[:pfam])

    outfile = "$(pfam).$(pdb).$(chain).csv"

    output = DataFrame(
      threshold = vcat([ 1.0, 0.62 ], rand(Beta(2,2), 13)),
      beta      = vcat([ 0.0, 4.60 ], rand(Beta(1,4), 13) .* 50.00),
      lambda    = vcat([ 0.0, 0.05 ], rand(Beta(0.5), 13) .* 10.00),

      nclusters  =  fill!(Array(Int, 15),-1),

      AUC_MIp    =  fill!(Array(Float64, 15), NaN),
      AUC_ZMIp   =  fill!(Array(Float64, 15), NaN),
      AUC_BLMIp  =  fill!(Array(Float64, 15), NaN),
      AUC_ZBLMIp =  fill!(Array(Float64, 15), NaN)
      )

    try
      msa, truecontacts, notcontacts, coverage = setup!(pfam, uniprot, pdb, chain, pdbdir, siftsdir)
      if (sum(truecontacts) != 0) && (sum(notcontacts) != 0)
        for row in eachrow(output)
          clustering = hobohmI(msa, row[:threshold])
          row[:nclusters] = getnclusters(clustering)
          zmip, mip = _buslje09(msa, row[:lambda], clustering)
          zblmip, blmip = _blmi(msa, row[:beta],   clustering)
          row[:AUC_MIp]    = AUC(mip,    truecontacts, notcontacts)
          row[:AUC_ZMIp]   = AUC(zmip,   truecontacts, notcontacts)
          row[:AUC_BLMIp]  = AUC(blmip,  truecontacts, notcontacts)
          row[:AUC_ZBLMIp] = AUC(zblmip, truecontacts, notcontacts)
        end
      else
        println("CMAP ERROR: $pfam $pdb $chain $uniprot truecontacts=$(sum(truecontacts)) notcontacts=$(sum(notcontacts))")
      end
      output[:coverage] = coverage
      output[:contacts]    = sum(truecontacts)
      output[:notcontacts] = sum(notcontacts)
    catch err
      println("ERROR FOR: ", name, " ",err)
    end

    output[:meanpid] = meanpercentidentity(msa)
    output[:ncol] = ncolumns(msa)
    output[:nseq] = nsequences(msa)

    output[:pfam]    = pfam
    output[:uniprot] = uniprot
    output[:pdb]     = pdb
    output[:chain]   = chain

    writetable(outfile, output)

    nothing
  end

end

pmap(optim_row, eachrow(data[1:3,:]))

