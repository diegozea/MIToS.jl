"""
`AUC(scores::PairwiseListMatrix, msacontacts::PairwiseListMatrix)`

Returns the Area Under a ROC (Receiver Operating Characteristic) Curve (AUC) of the
`scores` for `msacontact` prediction. `score` and `msacontact` lists are vinculated
(inner join) by their labels (i.e. column number in the file). `msacontact` should have
1.0 for true contacts and 0.0 for not contacts (NaN or other numbers for missing
values). You need to do `using ROCAnalysis` before using this function.
"""
function ROCAnalysis.AUC(scores::NamedArray{L,2,PairwiseListMatrix{L,false,VL},NL},
    msacontacts::NamedArray{R,2,PairwiseListMatrix{L,false,VR},NR}) where {L <: AbstractFloat,
        R <: AbstractFloat,VL,VR,NL,NR}
    sco, con = join(scores, msacontacts, kind=:inner)
    true_contacts, false_contacts = getcontactmasks(con)
    scores_list = getlist(getarray(sco))
    ROCAnalysis.AUC(
        ROCAnalysis.roc(
            scores_list[true_contacts  .& .!isnan.(scores_list)],
            scores_list[false_contacts .& .!isnan.(scores_list)]))
end

export AUC
