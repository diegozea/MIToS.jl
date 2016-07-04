const _residue_labels = map(string, reverse!(res"ARNDCQEGHILKMFPSTWYV-"))

@recipe function plot(msa::AbstractMatrix{Residue})
    seriestype --> :heatmap
    yflip --> true
    grid --> false
    foreground_color_border --> nothing
    foreground_color_axis --> nothing
    linewidth --> 0
    zdiscrete_values --> _residue_labels
    nseq = nsequences(msa)
    if nseq > 20
        step = div(nseq,20)
        yticks --> (1:step:nseq,names(msa))
        html_output_format :=  :png
    end
    1:ncolumns(msa), names(msa), map(string, msa)
end
