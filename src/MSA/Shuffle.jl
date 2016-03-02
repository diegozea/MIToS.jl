"""
Shuffles the residues in each column, keeping fixed the gap positions
"""
function shuffle_residues_columnwise!(aln::Matrix{Residue})
    nseq, nres = size(aln)
    for i in 1:nres
        @inbounds for j in 1:nseq
            a = aln[j,i]
            if a != GAP
                k = rand(1:nseq)
                b = aln[k,i]
                while b == GAP
                    k = rand(1:nseq)
                    b = aln[k,i]
                end
                aln[k,i] = a
                aln[j,i] = b
            end
        end
    end
    aln
end

# 0.00084 seconds faster than an implemetation similar to shuffle_residues_columnwise (PF00085)
"""
Shuffles the residues in each sequence, keeping fixed the gap positions
"""
function shuffle_residues_sequencewise!(aln::Matrix{Residue})
    taln = transpose(aln)
    shuffle_residues_columnwise!(taln)
    transpose!(aln, taln)
end

"""
Shuffles the residues in each sequence
"""
function shuffle_sequencewise!(aln::Matrix{Residue})
    nseq, nres = size(aln)
    for i in 1:nseq
        @inbounds for j in 1:nres
            k = rand(1:nres)
            aln[i,k], aln[i,j] = aln[i,j], aln[i,k]
        end
    end
    aln
end

"""
Shuffles the residues in each column
"""
function shuffle_columnwise!(aln::Matrix{Residue})
    nseq, nres = size(aln)
    for i in 1:nres
        @inbounds for j in 1:nseq
            k = rand(1:nseq)
            aln[k,i], aln[j,i] = aln[j,i], aln[k,i]
        end
    end
    aln
end

for fun in [ :shuffle_columnwise!, :shuffle_sequencewise!, :shuffle_residues_sequencewise!, :shuffle_residues_columnwise! ]
    @eval $(fun)(aln::AbstractMultipleSequenceAlignment) = $(fun)(aln.msa)
end
