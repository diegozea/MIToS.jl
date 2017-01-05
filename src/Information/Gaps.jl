# Pairwise Gap Percentage
# =======================

function gap_intersection_percentage{T}(nxy::Counts{T,2,GappedAlphabet})
    T(100.0) * gettablearray(nxy)[21,21] / gettotal(nxy)
end

function gap_union_percentage{T}(nxy::Counts{T,2,GappedAlphabet})
    marginals = getmarginalsarray(nxy)
    T(100.0) * (marginals[21,1] + marginals[21,2] - gettablearray(nxy)[21,21]) / gettotal(nxy)
end


# MIToS Pairwise Gap Percentage
# =============================

# """
# This function takes a MSA or a file and a `Format` as first arguments.
# Calculates the percentage of gaps on columns pairs (union and intersection) using sequence clustering (Hobohm I).
#
# Argument, type, default value and descriptions:
#
#   - clustering  Bool      true    Sequence clustering (Hobohm I)
#   - threshold             62      Percent identity threshold for sequence clustering (Hobohm I)
#
# This function returns:
#
#   - pairwise gap percentage (union)
#   - pairwise gap percentage (intersection)
# """
function pairwisegapfraction(aln::AbstractMatrix{Residue}; clustering::Bool=true, threshold=62)
    clusters = clustering ? hobohmI(aln, threshold) : NoClustering()
    table = Counts(ContingencyTable(Float64,Val{2},GappedAlphabet()))
    gu = mapcolpairfreq!(gap_union_percentage, aln, table, Val{true}, weights=clusters)
    gi = mapcolpairfreq!(gap_intersection_percentage, aln, table, Val{true}, weights=clusters)
    gu, gi
end

function pairwisegapfraction{T <: Format}(filename::String, format::Type{T}; kargs...)
    aln = read(filename, T, AnnotatedMultipleSequenceAlignment, generatemapping=true)
    pairwisegapfraction(aln; kargs...)
end
