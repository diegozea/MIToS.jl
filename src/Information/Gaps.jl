# Pairwise Gap Percentage
# =======================

"It calculates the gap intersection as percentage from a table of `Counts`."
function gap_intersection_percentage(nxy::Counts{T,2,GappedAlphabet}) where T
    T(100.0) * gettablearray(nxy)[21,21] / gettotal(nxy)
end

"It calculates the gap union as percentage from a table of `Counts`."
function gap_union_percentage(nxy::Counts{T,2,GappedAlphabet}) where T
    marginals = getmarginalsarray(nxy)
    T(100.0) * (marginals[21,1] + marginals[21,2] - gettablearray(nxy)[21,21]) / gettotal(nxy)
end


# MIToS Pairwise Gap Percentage
# =============================

"""
It takes a MSA or a file and a `FileFormat` as first arguments. It calculates the percentage
of gaps on columns pairs (union and intersection) using sequence clustering (Hobohm I).

Argument, type, default value and descriptions:

```
    - clustering  Bool      true    Sequence clustering (Hobohm I)
    - threshold             62      Percent identity threshold for sequence clustering (Hobohm I)
```

This function returns:

```
    - pairwise gap union as percentage
    - pairwise gap intersection as percentage
```  
"""
function pairwisegapfraction(aln::AbstractMatrix{Residue}; clustering::Bool=true, threshold=62)
    clusters = clustering ? hobohmI(aln, threshold) : NoClustering()
    table = Counts(ContingencyTable(Float64,Val{2},GappedAlphabet()))
    gu = mapcolpairfreq!(gap_union_percentage, aln, table, Val{true}, weights=clusters)
    gi = mapcolpairfreq!(gap_intersection_percentage, aln, table, Val{true}, weights=clusters)
    gu, gi
end

function pairwisegapfraction(filename::String, format::Type{T}; kargs...) where T <: FileFormat
    aln = read(filename, T, AnnotatedMultipleSequenceAlignment, generatemapping=true)
    pairwisegapfraction(aln; kargs...)
end
