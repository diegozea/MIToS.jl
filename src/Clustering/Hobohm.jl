# Hobohm I
# ========

"""
Fill `cluster` and `clustersize` matrices.
These matrices are assumed to be empty (only zeroes) and their length is assumed to be equal to the number
of sequences in the alignment (`aln`).
`threshold` is the minimum identity value between two sequences to be in the same cluster.
"""
function _fill_hobohmI!(cluster::Vector{Int}, clustersize::Vector{Int}, aln::Vector{Vector{Residue}}, threshold::Float64)
    cluster_id = 0
    nseq = length(aln)
    for i in 1:(nseq-1)
        if cluster[i] == 0
            cluster_id += 1
            cluster[i] = cluster_id
            clustersize[cluster_id] += 1
            ref_seq = aln[i]
            for j in (i+1):nseq
                if cluster[j] == 0 && percentidentity(ref_seq, aln[j], threshold)
                    cluster[j] = cluster_id
                    clustersize[cluster_id] += 1
                end
            end
        end
    end
    if cluster[nseq] == 0
        cluster_id += 1
        cluster[nseq] = cluster_id
        clustersize[cluster_id] += 1
    end
    resize!(clustersize, cluster_id)
end

"""
Calculates the weight of each sequence in a cluster.
The weight is equal to one divided by the number of sequences in the cluster.
"""
function _get_sequence_weight(clustersize, cluster)
    nseq = length(cluster)
    sequence_weight = Array(Float64, nseq)
    for i in 1:nseq
        @inbounds sequence_weight[i] = 1.0 / clustersize[cluster[i]]
    end
    sequence_weight
end

"""
Cluster sequences using Hobohm method 1.

For more information read:

1. Protein Sci. 1992 Mar;1(3):409-17.

Selection of representative protein data sets.

Hobohm U(1), Scharf M, Schneider R, Sander C.

Author information:
(1)European Molecular Biology Laboratory, Heidelberg, Germany.

The Protein Data Bank currently contains about 600 data sets of three-dimensional
protein coordinates determined by X-ray crystallography or NMR. There is
considerable redundancy in the data base, as many protein pairs are identical or
very similar in sequence. However, statistical analyses of protein
sequence-structure relations require nonredundant data. We have developed two
algorithms to extract from the data base representative sets of protein chains
with maximum coverage and minimum redundancy. The first algorithm focuses on
optimizing a particular property of the selected proteins and works by successive
selection of proteins from an ordered list and exclusion of all neighbors of each
selected protein. The other algorithm aims at maximizing the size of the selected
set and works by successive thinning out of clusters of similar proteins. Both
algorithms are generally applicable to other data bases in which criteria of
similarity can be defined and relate to problems in graph theory. The largest
nonredundant set extracted from the current release of the Protein Data Bank has
155 protein chains. In this set, no two proteins have sequence similarity higher
than a certain cutoff (30% identical residues for aligned subsequences longer
than 80 residues), yet all structurally unique protein families are represented.
Periodically updated lists of representative data sets are available by
electronic mail from the file server netserv@embl-heidelberg.de. The selection
may be useful in statistical approaches to protein folding as well as in the
analysis and documentation of the known spectrum of three-dimensional protein
structures.

PMCID: PMC2142204
PMID: 1304348  [PubMed - indexed for MEDLINE]
"""
function hobohmI(msa::AbstractMatrix{Residue}, threshold::Float64)
    aln = getresiduesequences(msa)
    nseq = length(aln)
    cluster = zeros(Int,nseq)
    clustersize = zeros(Int,nseq)
    _fill_hobohmI!(cluster, clustersize, aln, threshold)
    Clusters(clustersize, cluster, _get_sequence_weight(clustersize, cluster))
end
