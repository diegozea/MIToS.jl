# # Linking structural and evolutionary information
#
#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__cookbook/notebooks/02_Linking_structural_and_evolutionary_information.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__cookbook/notebooks/02_Linking_structural_and_evolutionary_information.ipynb)
#
#
# ## Problem description
#
# It is a very common task to map sequence to structure residue number. For
# example, to link structural information coming from *PDB* and evolutionary
# information calculated from multiple sequence alignments.
#
# The naive way of mapping sequence and structure is to perform global pairwise
# alignment between the sequence and the *PDB* sequence (using the residues in
# *ATOM*). The problem with this approach is that the sequences can have
# missing regions and standard pairwise alignment algorithms often yield
# incorrect assignations around those regions
# [(Velankar et.al. 2013)](https://doi.org/10.1093/nar/gks1258). This is
# particularly important when aligning *PDB* sequences, that can have missing
# residues, and sequences coming from multiple sequence alignments, that can be
# incomplete or have unaligned regions (e.g. insert states).
#
# The [SIFTS](https://www.ebi.ac.uk/pdbe/docs/sifts/index.html)
# (Structure Integration with Function, Taxonomy and Sequences) database solves
# this problem and provides residue level mapping between PDB and other
# databases (e.g. *UniProt* and *Pfam*).
#
# The `SIFTS` module of MIToS has functions to access this residue level mapping
# between *PDB* and other databases. Also, MIToS keeps track of the residue
# number of each residue in a multiple sequence alignment (MSA) using its annotations.
# Both things together, allow the correct mapping of sequence and structure
# without performing error-prone pairwise alignments.
#
# Particular solutions depend on problem details, here we show some common ways
# to use MIToS and *SIFTS* to map evolutionary information calculated in an MSA
# (e.g. entropy) with structural information (e.g. B-factors).
#
# ## PDB and Pfam alignment mapping
#
# This is the easiest problem to solve with the
# [MIToS `Pfam` module](@ref Module-Pfam) because *SIFTS* already has a residue
# level mapping between *PDB* and *Pfam*.
#
# For this example, we are going to map the columns in the multiple sequence
# alignment of the *PF09645* Pfam family and the residues in the chain *A* from
# the *2VQC* *PDB* file. The needed files are available in the MIToS test suite:

using MIToS
pdb_file   = abspath(pathof(MIToS), "..", "..", "test", "data", "2VQC.pdb")
pfam_file  = abspath(pathof(MIToS), "..", "..", "test", "data", "PF09645_full.stockholm")
sifts_file = abspath(pathof(MIToS), "..", "..", "test", "data", "2vqc.xml.gz")
#md nothing # hide

# You can also use `downloadpdb` from `MIToS.PDB`, `downloadpfam` from
# `MIToS.Pfam` and `downloadsifts` from `MIToS.SIFTS` to get the corresponding
# files from those databases.
#
# It is important to read the Pfam MSA file using `generatemapping=true` and
# `useidcoordinates=true` because that allows keeping track of the residue
# number using the MSA annotations.

using MIToS.Pfam
msa = read(pfam_file, Stockholm, generatemapping=true, useidcoordinates=true)

# First, we need to know what is the sequence in the MSA that correspond to the
# PDB we want to link. Luckily, Pfam Stockholm files store the mapping between
# sequences and PDB chains. You can access that mapping using the `getseq2pdb`
# function from `MIToS.Pfam`

seq2pdbs = getseq2pdb(msa)

# The returned dictionary gives you all the PDB chains associated with a
# determined sequence in the MSA. But, in this case, we want to go in the other
# direction to find all the sequences associated with a determined *PDB* chain.
# We are going to use a list comprehension because it is possible for a single
# chain to be associated with more than one sequence in the *Pfam* MSA (e.g.
# domain repeats).

pdb_code  = "2VQC"
pdb_chain = "A"
seq_ids = [ seq for (seq, pdbs) in seq2pdbs if (pdb_code, pdb_chain) in pdbs ]

# In this example, we are going to use the only sequence we found for the *A*
# of *2VQC*.

seq_id = seq_ids[1]

# Finally, we can use the `msacolumn2pdbresidue` function from the Pfam module
# to get a dictionary from the MSA column index to the *PDB* residue number:

pfam_id = "PF09645"
msacol2pdbres = msacolumn2pdbresidue(msa, seq_id, pdb_code, pdb_chain, pfam_id, sifts_file)

# This dictionary has the mapping between MSA column and PDB residue that allows
# the mapping between evolutionary and structural information. For example, to
# measure the correlation between entropy (related to residue variation in an
# MSA column) and the mean B factor of the residue:

using MIToS.Information
Hx = mapcolfreq!(entropy,
				 msa,
				 Counts(ContingencyTable(Int, Val{1}, UngappedAlphabet())))

# To get quick access to each PDB residue based on its residue number, we can
# read the PDB file into a dictionary using the `read` and `residuesdict`
# functions from the MIToS `PDB` module:

using MIToS.PDB
res_dict = residuesdict(read(pdb_file, PDBFile, occupancyfilter=true), "1", "A") # model 1 chain A

# Then, we can iterate the mapping dictionary to link the MSA and PDB based
# values:

using Statistics

x = Float64[]
y = Float64[]

for (col_index, res_number) in msacol2pdbres
	if res_number != "" # i.e. MSA column has an associated PDB residue
		push!(x, Hx[col_index])
		push!(y, mean(parse(Float64, atom.B) for atom in res_dict[res_number].atoms))
	end
end

cor(x, y)

# ##  Unknown sequence coordinates
#
# While *Pfam* alignments have the start and end of the aligned region indicated
# in the sequence name, other multiple sequence alignments don't give any hint
# about that. In those cases, we should use pairwise alignments. However,
# instead of aligning the sequence coming from the MSA and the *PDB* sequence,
# we can align the MSA sequence to the *UniProt* sequence to reduce the
# possibility of mapping errors. Once we have the mapping of the MSA sequence
# to the *UniProt* sequence, we can use *SIFTS* to map the *PDB* sequence to
# the MSA sequence using the *UniProt* numeration.
#
# For this example, we are going to use the following files included in MIToS
# documentation:

using MIToS
pdb_file   = abspath(pathof(MIToS), "..", "..", "docs", "data", "1dur.pdb")
msa_file   = abspath(pathof(MIToS), "..", "..", "docs", "data", "blast_alignment.fa")
sifts_file = abspath(pathof(MIToS), "..", "..", "docs", "data", "1dur.xml.gz")
uniprot_file = abspath(pathof(MIToS), "..", "..", "docs", "data", "P00193.fasta")
#md nothing # hide

# First, we are going to read the MSA file. In this case, we can not use
# `useidcoordinates=true` because the sequence names don't have the sequence
# coordinates in the Pfam format. However, we are going to use
# `generatemapping=true` to get the default mapping for each sequence in the
# alignment (from `1` to the length of the aligned region):

using MIToS.MSA
msa = read(msa_file, FASTA, generatemapping=true)

# After that, we get the first sequence of the MSA, the one we know that
# corresponds to the PDB of interest. We need the sequence as a `String`
# without gaps (unaligned), so we use the `MIToS.MSA` `stringsequence` function
# together with `replace`:

msa_seq = replace(stringsequence(msa, 1), '-' => "")

# Also, we are going to read the *UniProt* sequence. You can easily download the
# sequence from UniProt by doing:
# ```julia
# using MIToS.Utils
# download_file("https://www.uniprot.org/uniprot/P00193.fasta", "P00193.fasta")
# ```
# To read the FASTA file we are going to use the `FastaIO` package:

using FastaIO
uniprot_sequences = readfasta(uniprot_file)

# And get the unique sequence:

uniprot_seq = uniprot_sequences[1][2]

# We can perform a pairwise sequence alignment between both sequences by using
# the [`BioAlignments` package](https://github.com/BioJulia/BioAlignments.jl)
# from the *BioJulia* suite. In this case, we use a semi-global alignment
# (no start/end gap penalty) because we know that the MSA sequence is a region
# of the *UniProt* sequence.

using BioAlignments
costmodel = AffineGapScoreModel(BLOSUM62, gap_open=-10, gap_extend=-1)
aln = pairalign(SemiGlobalAlignment(), msa_seq, uniprot_seq, costmodel)

# Then, we only need to iterate the alignment to designate the positions and
# store the equivalences in a dictionary:

function seq2refnumber(aln)
    seq_pos = 0
	ref_pos = 0
	last_seq_pos = 0
	seq2ref = Dict{Int, Int}()
    for (seq_res, ref_res) in alignment(aln)
        if seq_res != '-'
            seq_pos += 1
		end
        if ref_res != '-'
            ref_pos += 1
		end
		if seq_pos != last_seq_pos
			seq2ref[seq_pos] = ref_pos
			last_seq_pos = seq_pos
    	end
	end
    seq2ref
end

seqnum2uniprotnum = seq2refnumber(aln)

# Then, we can use `getsequencemapping` to go from MSA column number to
# *UniProt* residue, and `siftsmapping` to go from *UniProt* to *PDB*:

seqmap = getsequencemapping(msa, 1)

#-

colnum2uniprotnum = Dict{Int, Int}()
for (colnum, seqnum) in enumerate(seqmap)
	if seqnum != 0 # getsequencemapping returns 0 where there is a gap
		colnum2uniprotnum[colnum] = seqnum2uniprotnum[seqnum]
	end
end
colnum2uniprotnum

#-

using MIToS.SIFTS

uniprotnum2pdbnum = siftsmapping(sifts_file,
    dbUniProt,
    "P00193",
    dbPDB,
    "1dur", # SIFTS stores PDB identifiers in lowercase
    chain="A",
    missings=false) # residues without coordinates aren't used in the mapping

# To finally get the dictionary from MSA column index to PDB residue number

colnum2pdbnum = Dict{Int, String}()
for (colnum, uniprotnum) in colnum2uniprotnum
	pdbresnum = get(uniprotnum2pdbnum, string(uniprotnum), "")
	if pdbresnum != ""
		colnum2pdbnum[colnum] = pdbresnum
	end
end

colnum2pdbnum
