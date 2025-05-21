const _residue_labels = map(string, reverse!(res"ARNDCQEGHILKMFPSTWYV-"))

"""
    plot(msa::AbstractMultipleSequenceAlignment)

Defines a plot recipe for `AbstractMultipleSequenceAlignment` types (and subtypes)
using the Plots.jl package. This allows direct plotting of MSA objects.

The recipe generates a heatmap representation of the MSA:
- **`seriestype`**: `:heatmap`
- **`yflip`**: `true` (sequences typically start from top)
- **`grid`**: `false`
- **`zdiscrete_values`**: Maps residues to string labels for hover/colorbar.
- **Ticks**: For MSAs with more than 20 sequences, y-axis ticks are subsampled for readability.
- **Output**: Residue data is converted to strings for heatmap categorization.

# Example
```julia
using Plots
using MIToS.MSA
# Assuming `msa_object` is an AbstractMultipleSequenceAlignment
# plot(msa_object)
```
"""
@recipe function plot(msa::AbstractMultipleSequenceAlignment)
    seriestype --> :heatmap
    yflip --> true
    grid --> false
    foreground_color_border --> nothing
    foreground_color_axis --> nothing
    linewidth --> 0
    zdiscrete_values --> _residue_labels
    nseq = nsequences(msa)
    names = sequencenames(msa)
    if nseq > 20
        step = div(nseq, 20)
        yticks --> (1:step:nseq, names)
        html_output_format := :png
    end
    1:ncolumns(msa), names, map(string, getresidues(msa))
end
